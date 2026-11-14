#!/usr/bin/env python3
"""测试夹具生成器：自写 BGZF/BAM 写出端（与被测实现相互独立）。

生成：
- fixtures/valid_*.bam      合法输入，附 *.expected.json（生成器自算投影+vo）
- fixtures/invalid_*.bam    非法输入，期望整次拒绝
- fixtures/manifest.json    每个夹具的期望（ok / error+offset）

用法：python3 tests/make_fixtures.py
"""

import json
import os
import random
import struct
import zlib

HERE = os.path.dirname(os.path.abspath(__file__))
FIX = os.path.join(HERE, "fixtures")

BASES = "=ACMGRSVTWYHKDBN"
OPS = "MIDNSHP=XB"


# ---------------------------------------------------------------- BGZF 写出端

def bgzf_block(payload, *, flg=4, n_bc=1, extra_tail=b"", raw_extra=None,
               corrupt=None):
    """构造一个 BGZF 块（含 BC 子字段与正确 BSIZE/CRC32/ISIZE）。"""
    co = zlib.compressobj(6, zlib.DEFLATED, -15)
    comp = co.compress(payload) + co.flush()

    def build(bsize):
        if raw_extra is not None:
            extra = raw_extra
        else:
            extra = b""
            for _ in range(n_bc):
                extra += b"BC" + struct.pack("<HH", 2, bsize)
            extra += extra_tail
        hdr = struct.pack("<BBBBIBBH", 0x1F, 0x8B, 8, flg, 0, 0, 0xFF,
                          len(extra)) + extra
        return hdr + comp + struct.pack(
            "<II", zlib.crc32(payload) & 0xFFFFFFFF, len(payload) & 0xFFFFFFFF)

    blk = build(0)
    blk = build(len(blk) - 1)  # BSIZE = 块总字节数 - 1
    if corrupt:
        blk = corrupt(blk)
    return blk


EOF_BLOCK = bgzf_block(b"")  # 末块须空


def bgzf_stream(payload, chunk=None, empty_after=()):
    """把解压流切成若干 BGZF 块；empty_after 指定在哪些块之后插入空块。

    返回 (compressed, blocks)，blocks 为 [(ustart, usize, cstart)]（仅非空块）。
    """
    if chunk is None:
        parts = [payload] if payload or True else []
    else:
        parts = [payload[i:i + chunk] for i in range(0, len(payload), chunk)]
    blocks = []
    out = bytearray()
    ustart = 0
    for i, part in enumerate(parts):
        blk = bgzf_block(part)
        if part:
            blocks.append((ustart, len(part), len(out)))
            ustart += len(part)
        out += blk
        if i in empty_after:
            out += bgzf_block(b"")  # 内部空块
    out += EOF_BLOCK
    return bytes(out), blocks


def vo_of(blocks, u):
    for us, sz, cs in blocks:
        if us <= u < us + sz:
            return (cs << 16) | (u - us)
    raise ValueError("offset %d not in any block" % u)


# ----------------------------------------------------------------- BAM 写出端

def pack_seq(seq, pad=0):
    codes = [BASES.index(c) for c in seq]
    out = bytearray()
    for i in range(0, len(codes), 2):
        hi = codes[i]
        lo = codes[i + 1] if i + 1 < len(codes) else pad
        out.append((hi << 4) | lo)
    return bytes(out)


def cigar_bytes(ops):
    return b"".join(
        struct.pack("<I", (ln << 4) | (OPS.index(op) if isinstance(op, str) else op))
        for ln, op in ops)


def aux_bytes(fields):
    out = bytearray()
    for tag, ty, val in fields:
        out += tag.encode("ascii") + ty.encode("ascii")
        if ty == "A":
            out += val.encode("ascii")
        elif ty == "c":
            out += struct.pack("<b", val)
        elif ty == "C":
            out += struct.pack("<B", val)
        elif ty == "s":
            out += struct.pack("<h", val)
        elif ty == "S":
            out += struct.pack("<H", val)
        elif ty == "i":
            out += struct.pack("<i", val)
        elif ty == "I":
            out += struct.pack("<I", val)
        elif ty == "f":
            out += struct.pack("<f", val)
        elif ty in ("Z", "H"):
            out += val.encode("ascii") + b"\0"
        elif ty == "B":
            sub, vals = val
            fmt = {"c": "b", "C": "B", "s": "h", "S": "H",
                   "i": "i", "I": "I", "f": "f"}[sub]
            out += sub.encode("ascii") + struct.pack("<i", len(vals))
            out += struct.pack("<%d%s" % (len(vals), fmt), *vals)
        else:
            raise ValueError(ty)
    return bytes(out)


def reg2bin(beg, end):
    end -= 1
    if beg >> 14 == end >> 14:
        return 4681 + (beg >> 14)
    if beg >> 17 == end >> 17:
        return 585 + (beg >> 17)
    if beg >> 20 == end >> 20:
        return 73 + (beg >> 20)
    if beg >> 23 == end >> 23:
        return 9 + (beg >> 23)
    if beg >> 26 == end >> 26:
        return 1 + (beg >> 26)
    return 0


def record(ref_id, pos, name, flag, cigar, seq, qual, aux=b"", mapq=60,
           pad=0, name_nul=True, tail=b"", l_seq=None, next_ref=-1, next_pos=-1):
    rn = name.encode("utf-8") + (b"\0" if name_nul else b"")
    cig = cigar if isinstance(cigar, bytes) else cigar_bytes(cigar)
    sq = pack_seq(seq, pad) if isinstance(seq, str) else seq
    ls = len(seq) if l_seq is None else l_seq
    ql = bytes(qual)
    span = 0
    for ln, op in (cigar if not isinstance(cigar, bytes) else []):
        if (OPS.index(op) if isinstance(op, str) else op) in (0, 2, 3, 7, 8):
            span += ln
    bin_ = reg2bin(pos, pos + max(span, 1)) if ref_id >= 0 and pos >= 0 else 0
    core = struct.pack("<iiBBHHHiiii", ref_id, pos, len(rn), mapq, bin_,
                       len(cig) // 4, flag, ls, next_ref, next_pos, 0)
    body = core + rn + cig + sq + ql + aux + tail
    return struct.pack("<i", len(body)) + body


def header_stream(refs, text=b""):
    out = b"BAM\x01" + struct.pack("<i", len(text)) + text
    out += struct.pack("<i", len(refs))
    for name, lref in refs:
        nb = name.encode("utf-8") + b"\0"
        out += struct.pack("<i", len(nb)) + nb + struct.pack("<i", lref)
    return out


def project(pos, cigar, seq):
    """生成器侧独立投影（与被测实现相互独立）。"""
    cols = []
    r, q = pos, 0
    for ln, op in cigar:
        o = OPS.index(op) if isinstance(op, str) else op
        if o in (0, 7, 8):
            for _ in range(ln):
                cols.append([r, q, seq[q]])
                r += 1
                q += 1
        elif o in (1, 4):
            q += ln
        elif o in (2, 3):
            r += ln
    return cols


# ------------------------------------------------------------- 合法夹具的内容

REFS = [("chr1", 1000), ("chr2", 50)]
TEXT = (b"@HD\tVN:1.6\tSO:unknown\n"
        b"@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:50\n")

# (name, ref_id, pos, flag, cigar, seq, qual, aux, pad)
MAIN_RECORDS = [
    # 任务书示例：pos=0、SEQ=ACG、CIGAR=1M1I1M → [[0,0,"A"],[1,2,"G"]]
    ("r001", 0, 0, 0, [(1, "M"), (1, "I"), (1, "M")], "ACG", [30, 30, 30], b"", 0),
    ("r002", 0, 10, 0, [(3, "M")], "ACG", [255, 255, 255], b"", 0),
    # 负链：不再反补，用存储 SEQ
    ("r003", 0, 20, 16, [(2, "S"), (4, "M")], "TTACGN", [20] * 6, b"", 0),
    ("r004", 0, 30, 0, [(2, "="), (1, "X"), (2, "D"), (1, "M")], "ACGT",
     [255] * 4, b"", 0),
    ("r005", 0, 40, 0, [(2, "M"), (3, "N"), (2, "M")], "ACGT", [31] * 4, b"", 0),
    ("r006", 0, 50, 0, [(2, "M"), (2, "P"), (2, "M")], "ACGT", [32] * 4, b"", 0),
    ("r007", 0, 60, 0, [(1, "H"), (2, "S"), (3, "M"), (2, "S"), (1, "H")],
     "TTACGAA", [33] * 7, b"", 0),
    # 未比对：跳过
    ("r008", -1, -1, 4, [], "ACG", [10, 10, 10], b"", 0),
    # 已比对但 SEQ 为空：跳过
    ("r009", 0, 5, 0, [], "", [], b"", 0),
    # 次比对 / 补充比对：仍收录
    ("r010", 1, 0, 0x100, [(1, "M")], "N", [40], b"", 0),
    ("r011", 1, 10, 0x800, [(2, "M")], "RY", [41, 41], b"", 0),
    # 重名保留
    ("r001", 1, 20, 0, [(1, "M")], "=", [93], b"", 0),
    # 全部辅助类型
    ("r012", 0, 100, 0, [(5, "M")], "ACGTN", [25] * 5,
     aux_bytes([
         ("AA", "A", "x"), ("BC", "c", -5), ("BD", "C", 200),
         ("BE", "s", -3000), ("BF", "S", 60000),
         ("BG", "i", -70000), ("BH", "I", 4000000000),
         ("BI", "f", 3.5), ("BJ", "Z", "hello"), ("BK", "H", "DEADBEEF"),
         ("BL", "B", ("c", [-1, 2])), ("BM", "B", ("C", [1, 255])),
         ("BN", "B", ("s", [-1000, 2])), ("BO", "B", ("S", [1, 65535])),
         ("BP", "B", ("i", [-100000, 2])), ("BQ", "B", ("I", [1, 4294967295])),
         ("BR", "B", ("f", [1.5, -2.5])),
     ]), 0),
    # 奇数长度 SEQ，末半字节填 0xF（应被忽略）
    ("r013", 0, 110, 0, [(5, "M")], "ACGTA", [30] * 5, b"", 0xF),
    ("r014", 0, 120, 0, [(2, "M")], "AC", [0, 93], b"", 0),
    # 末碱基对齐（跨度恰好到参考末端）
    ("r015", 1, 49, 0, [(1, "M")], "A", [30], b"", 0),
    # 带 FLAG 0x4 但写了 refID/pos：仍按未比对跳过
    ("r016", 0, 0, 4, [], "AC", [30, 30], b"", 0),
    # 末尾 D
    ("r017", 0, 130, 0, [(2, "M"), (1, "D")], "AC", [30, 30], b"", 0),
]


def main_records_stream(text=TEXT):
    stream = header_stream(REFS, text)
    expect_meta = []  # (ustart, name, ref, columns)
    for (name, rid, pos, flag, cig, seq, qual, aux, pad) in MAIN_RECORDS:
        u = len(stream)
        stream += record(rid, pos, name, flag, cig, seq, qual, aux=aux, pad=pad)
        if not (flag & 0x4) and len(seq) > 0:
            expect_meta.append((u, name, REFS[rid][0], project(pos, cig, seq)))
    return stream, expect_meta


def build_valid(name, stream, expect_meta, chunk=None, empty_after=()):
    comp, blocks = bgzf_stream(stream, chunk=chunk, empty_after=empty_after)
    expected = [{"name": nm, "ref": rf, "vo": str(vo_of(blocks, u)),
                 "columns": cols} for u, nm, rf, cols in expect_meta]
    write(name + ".bam", comp)
    write(name + ".expected.json",
          json.dumps(expected, separators=(",", ":")).encode("utf-8"))
    return {"status": "ok", "bam": name + ".bam",
            "expected": name + ".expected.json"}


# ------------------------------------------------------------------ 随机合法

def random_records_stream(seed=20261114, count=60):
    rng = random.Random(seed)
    refs = [("chrR", 100000)]
    text = b"@HD\tVN:1.6\n@SQ\tSN:chrR\tLN:100000\n"
    stream = header_stream(refs, text)
    meta = []
    for i in range(count):
        ops = []
        lead_h = rng.randint(0, 2) if rng.random() < 0.3 else 0
        lead_s = rng.randint(1, 3) if rng.random() < 0.5 else 0
        trail_s = rng.randint(1, 3) if rng.random() < 0.5 else 0
        trail_h = rng.randint(0, 2) if rng.random() < 0.3 else 0
        if lead_h:
            ops.append((lead_h, "H"))
        if lead_s:
            ops.append((lead_s, "S"))
        consume = lead_s + trail_s
        span = 0
        for _ in range(rng.randint(1, 6)):
            o = rng.choice("MIDN=XP")
            ln = rng.randint(1, 8)
            ops.append((ln, o))
            if o in "MI=X":
                consume += ln
            if o in "MDN=X":
                span += ln
        if trail_s:
            ops.append((trail_s, "S"))
        if trail_h:
            ops.append((trail_h, "H"))
        if consume == 0:
            ops.append((1, "M"))
            consume += 1
            span += 1
        seq = "".join(rng.choice("=ACMGRSVTWYHKDBN") for _ in range(consume))
        pos = rng.randint(0, 100000 - span)
        flag = rng.choice([0, 0, 16, 0x100, 0x800, 16 | 0x100])
        qual_mode = rng.random()
        if qual_mode < 0.4:
            qual = [255] * consume
        elif qual_mode < 0.8:
            qual = [rng.randint(0, 93) for _ in range(consume)]
        else:
            qual = [30] * consume
        aux = b""
        if rng.random() < 0.5:
            aux = aux_bytes([("NM", "i", rng.randint(0, 5)),
                             ("RG", "Z", "grp%d" % rng.randint(0, 3))])
        name = "rnd%03d" % i
        u = len(stream)
        stream += record(0, pos, name, flag, ops, seq, qual, aux=aux)
        meta.append((u, name, "chrR", project(pos, ops, seq)))
    # 末尾追加未比对记录
    u = len(stream)
    stream += record(-1, -1, "rnd_unmapped", 4, [], "ACGT", [255] * 4)
    return stream, meta


# ------------------------------------------------------------- CG 占位 CIGAR

def cg_stream():
    refs = [("chrL", 200000)]
    text = b"@HD\tVN:1.6\n@SQ\tSN:chrL\tLN:200000\n"
    stream = header_stream(refs, text)
    # 真实 CIGAR：(1M1D) x 32768，共 65536 个操作 > 65535 → 占位 kSmN
    real_ops = [(1, "M"), (1, "D")] * 32768
    seq = "ACGT" * 8192  # l_seq = 32768
    cg_vals = []
    for ln, op in real_ops:
        cg_vals.append((ln << 4) | OPS.index(op))
    aux = aux_bytes([("CG", "B", ("I", cg_vals))])
    placeholder = [(len(seq), "S"), (65535, "N")]
    u = len(stream)
    stream += record(0, 0, "cg1", 0, placeholder, seq, [255] * len(seq), aux=aux)
    meta = [(u, "cg1", "chrL", project(0, real_ops, seq))]
    return stream, meta


# ------------------------------------------------------------------- 写文件

def write(name, data):
    with open(os.path.join(FIX, name), "wb") as fh:
        fh.write(data)


def corrupt_byte(idx):
    def _c(blk):
        b = bytearray(blk)
        b[idx] ^= 0xFF
        return bytes(b)
    return _c


def set_bsize(value):
    def _c(blk):
        b = bytearray(blk)
        b[16:18] = value.to_bytes(2, "little")
        return bytes(b)
    return _c


def main():
    os.makedirs(FIX, exist_ok=True)
    manifest = {}

    # ---------------- 合法夹具 ----------------
    stream, meta = main_records_stream()
    manifest["valid_main"] = build_valid("valid_main", stream, meta)
    # 同一内容切小块（记录跨块）并插入内部空块
    manifest["valid_multiblock"] = build_valid(
        "valid_multiblock", stream, meta, chunk=97, empty_after=(0, 3))
    # 随机合法记录
    stream, meta = random_records_stream()
    manifest["valid_random"] = build_valid("valid_random", stream, meta)
    # CG 占位 CIGAR（真实 CIGAR 65536 个操作，记录跨多块）
    stream, meta = cg_stream()
    manifest["valid_cg"] = build_valid("valid_cg", stream, meta, chunk=60000)
    # 无记录 → []
    stream = header_stream(REFS, TEXT)
    manifest["valid_empty"] = build_valid("valid_empty", stream, [])
    # 全部未比对 → []
    stream = header_stream(REFS, TEXT)
    stream += record(-1, -1, "u1", 4, [], "ACG", [255, 255, 255])
    stream += record(0, 0, "u2", 4, [], "AC", [30, 30])
    manifest["valid_all_unmapped"] = build_valid("valid_all_unmapped", stream, [])
    # 无文本头
    stream, meta = main_records_stream(text=b"")
    manifest["valid_no_text"] = build_valid("valid_no_text", stream, meta)

    # ---------------- 非法：BGZF 层（报压缩块零基起点） ----------------
    hdr = header_stream(REFS, TEXT)
    good = bgzf_block(hdr + record(0, 0, "x", 0, [(1, "M")], "A", [30]))
    good2 = bgzf_block(hdr)

    def invalid_bgzf(name, data, offset):
        write(name + ".bam", data)
        manifest[name] = {"status": "error", "bam": name + ".bam",
                          "offset": offset}

    invalid_bgzf("bgzf_empty_file", b"", 0)
    invalid_bgzf("bgzf_bad_magic", corrupt_byte(0)(good) + EOF_BLOCK, 0)
    invalid_bgzf("bgzf_bad_cm", corrupt_byte(2)(good) + EOF_BLOCK, 0)
    invalid_bgzf("bgzf_flg0", bgzf_block(hdr, flg=0) + EOF_BLOCK, 0)
    invalid_bgzf("bgzf_flg6", bgzf_block(hdr, flg=6) + EOF_BLOCK, 0)
    invalid_bgzf("bgzf_two_bc", bgzf_block(hdr, n_bc=2) + EOF_BLOCK, 0)
    invalid_bgzf("bgzf_no_bc",
                 bgzf_block(hdr, n_bc=0, extra_tail=b"XY\x02\x00ab") + EOF_BLOCK,
                 0)
    invalid_bgzf("bgzf_bc_slen4",
                 bgzf_block(hdr, raw_extra=b"BC" + struct.pack("<H", 4)
                            + struct.pack("<I", 0)) + EOF_BLOCK, 0)
    invalid_bgzf("bgzf_bad_crc", corrupt_byte(-8)(good) + EOF_BLOCK, 0)
    invalid_bgzf("bgzf_bad_isize", corrupt_byte(-4)(good) + EOF_BLOCK, 0)
    invalid_bgzf("bgzf_bad_deflate", corrupt_byte(20)(good) + EOF_BLOCK, 0)
    invalid_bgzf("bgzf_bsize_small", set_bsize(10)(good) + EOF_BLOCK, 0)
    invalid_bgzf("bgzf_bsize_overrun", set_bsize(60000)(good) + EOF_BLOCK, 0)
    invalid_bgzf("bgzf_isize_over_64k",
                 bgzf_block(b"\x00" * 70000) + EOF_BLOCK, 0)
    # 末块被截断：报该块起点
    invalid_bgzf("bgzf_truncated_tail", good + EOF_BLOCK[:-5], len(good))
    # 末块非空：两个数据块无 EOF，报第二块起点
    invalid_bgzf("bgzf_no_eof", good2 + good2, len(good2))
    # EOF 后垃圾：报垃圾起点
    invalid_bgzf("bgzf_garbage_after_eof",
                 good2 + EOF_BLOCK + b"\xde\xad\xbe\xef",
                 len(good2) + len(EOF_BLOCK))
    # 第二块 CRC 错：报第二块起点
    invalid_bgzf("bgzf_second_block_crc",
                 good2 + corrupt_byte(-8)(good2) + EOF_BLOCK, len(good2))

    # ---------------- 非法：头部（报 0） ----------------
    def invalid_hdr(name, stream):
        comp = bgzf_block(stream) + EOF_BLOCK
        write(name + ".bam", comp)
        manifest[name] = {"status": "error", "bam": name + ".bam", "offset": 0}

    invalid_hdr("hdr_bad_magic", b"BAM\x02" + hdr[4:])
    invalid_hdr("hdr_short", b"BAM")
    invalid_hdr("hdr_ltext_negative", b"BAM\x01" + struct.pack("<i", -1))
    invalid_hdr("hdr_ltext_overrun",
                b"BAM\x01" + struct.pack("<i", 1000) + b"short")
    invalid_hdr("hdr_nref_negative",
                b"BAM\x01" + struct.pack("<ii", 0, -1))
    invalid_hdr("hdr_nref_overrun",
                b"BAM\x01" + struct.pack("<ii", 0, 2)
                + struct.pack("<i", 5) + b"chr1\0" + struct.pack("<i", 100))
    invalid_hdr("hdr_ref_dup", header_stream([("chr1", 100), ("chr1", 200)]))
    invalid_hdr("hdr_ref_len_zero", header_stream([("chr1", 0)]))
    invalid_hdr("hdr_ref_len_negative", header_stream([("chr1", -5)]))
    invalid_hdr("hdr_ref_name_no_nul",
                b"BAM\x01" + struct.pack("<ii", 0, 1)
                + struct.pack("<i", 4) + b"chr1" + struct.pack("<i", 100))
    invalid_hdr("hdr_ref_lname_zero",
                b"BAM\x01" + struct.pack("<ii", 0, 1)
                + struct.pack("<i", 0) + struct.pack("<i", 100))
    invalid_hdr("hdr_ref_truncated",
                b"BAM\x01" + struct.pack("<ii", 0, 1)
                + struct.pack("<i", 10) + b"chr")

    # ---------------- 非法：记录（报 vo） ----------------
    def invalid_rec(name, recs, chunk=None):
        stream = hdr
        starts = []
        for r in recs:
            starts.append(len(stream))
            stream += r
        comp, blocks = bgzf_stream(stream, chunk=chunk)
        write(name + ".bam", comp)
        manifest[name] = {"status": "error", "bam": name + ".bam",
                          "offset": vo_of(blocks, starts[-1])}

    invalid_rec("rec_bad_op", [record(0, 0, "x", 0, [(1, 12)], "A", [30])])
    invalid_rec("rec_op_b", [record(0, 0, "x", 0, [(1, 9)], "A", [30])])
    invalid_rec("rec_zero_len", [record(0, 0, "x", 0, [(0, "M"), (1, "M")],
                                        "A", [30])])
    invalid_rec("rec_consumption",
                [record(0, 0, "x", 0, [(3, "M")], "ACGT", [30] * 4)])
    invalid_rec("rec_span", [record(1, 49, "x", 0, [(2, "M")], "AC", [30, 30])])
    invalid_rec("rec_pos_beyond",
                [record(1, 60, "x", 0, [(1, "M")], "A", [30])])
    invalid_rec("rec_soft_middle",
                [record(0, 0, "x", 0, [(2, "M"), (2, "S"), (2, "M")],
                        "ACACGT", [30] * 6)])
    invalid_rec("rec_hard_middle",
                [record(0, 0, "x", 0, [(2, "M"), (2, "H"), (2, "M")],
                        "ACGT", [30] * 4)])
    invalid_rec("rec_h_before_s_end",
                [record(0, 0, "x", 0, [(2, "M"), (1, "H"), (1, "S")],
                        "ACG", [30] * 3)])
    invalid_rec("rec_qual_mixed",
                [record(0, 0, "x", 0, [(3, "M")], "ACG", [255, 30, 255])])
    invalid_rec("rec_qual_94",
                [record(0, 0, "x", 0, [(2, "M")], "AC", [94, 30])])
    invalid_rec("rec_dup_aux",
                [record(0, 0, "x", 0, [(1, "M")], "A", [30],
                        aux=aux_bytes([("NM", "i", 1), ("NM", "i", 2)]))])
    invalid_rec("rec_bad_aux_type",
                [record(0, 0, "x", 0, [(1, "M")], "A", [30], aux=b"ZZQ\x00")])
    invalid_rec("rec_aux_truncated",
                [record(0, 0, "x", 0, [(1, "M")], "A", [30], aux=b"NMZabc")])
    invalid_rec("rec_aux_b_bad_sub",
                [record(0, 0, "x", 0, [(1, "M")], "A", [30],
                        aux=b"NMB" + b"x" + struct.pack("<i", 0))])
    invalid_rec("rec_aux_b_overrun",
                [record(0, 0, "x", 0, [(1, "M")], "A", [30],
                        aux=b"NMBi" + struct.pack("<ii", 10, 1))])
    invalid_rec("rec_aux_extra_byte",
                [record(0, 0, "x", 0, [(1, "M")], "A", [30], tail=b"\x00")])
    invalid_rec("rec_placeholder_no_cg",
                [record(0, 0, "x", 0, [(4, "S"), (65535, "N")], "ACGT",
                        [30] * 4)])
    invalid_rec("rec_placeholder_cg_z",
                [record(0, 0, "x", 0, [(4, "S"), (65535, "N")], "ACGT",
                        [30] * 4, aux=aux_bytes([("CG", "Z", "foo")]))])
    invalid_rec("rec_refid_range", [record(7, 0, "x", 0, [(1, "M")], "A", [30])])
    invalid_rec("rec_pos_negative",
                [record(0, -1, "x", 0, [(1, "M")], "A", [30])])
    invalid_rec("rec_name_no_nul",
                [record(0, 0, "x", 0, [(1, "M")], "A", [30], name_nul=False)])
    invalid_rec("rec_block_size_small", [struct.pack("<i", 20) + b"\x00" * 20])
    invalid_rec("rec_block_size_overrun",
                [struct.pack("<i", 1000000) + b"\x00" * 40])
    # 记录后残留 2 字节：报残留处 vo
    stream = hdr + record(0, 0, "ok", 0, [(1, "M")], "A", [30])
    stray = len(stream)
    stream += b"\xff\xfe"
    comp, blocks = bgzf_stream(stream)
    write("rec_leftover.bam", comp)
    manifest["rec_leftover"] = {"status": "error", "bam": "rec_leftover.bam",
                                "offset": vo_of(blocks, stray)}
    # 第二条记录出错且位于非首块：vo 含非零压缩起点
    good_rec = record(0, 0, "ok", 0, [(1, "M")], "A", [30])
    bad_rec = record(0, 0, "bad", 0, [(1, 15)], "A", [30])
    invalid_rec("rec_second_record_bad_multiblock", [good_rec, bad_rec],
                chunk=64)

    write("manifest.json", json.dumps(manifest, indent=1).encode("utf-8"))
    print("fixtures written to %s (%d entries)" % (FIX, len(manifest)))


if __name__ == "__main__":
    main()
