"""BAM 解析与参考投影（SAM/BAM 1.6 第 4 节）。自写实现，不依赖 BAM 读取器。

流程：先解析二进制头（文本头跳过；参考名唯一、长度正），再按文件序逐条
解析比对记录，仅收"已比对且 SEQ 非空"的记录（FLAG 0x4 或 l_seq==0 跳过），
对收录记录做全部语义校验并投影为 columns = [[r, q, b], ...]：

- 仅 M、=、X 产生项；r 为零基参考位置，q 为零基存储 SEQ 下标，b 为该字母；
- I、S 仅推进 q；D、N 仅推进 r；H、P 均不推进；负链不再反补；
- CIGAR 长度须正，操作码仅收 0..8（MIDNSHP=X，B 及其余非法）；
- 消费 SEQ 总长须等于 l_seq；参考跨度（M/D/N/=/X）不得越界；
- 剪切位置遵循规范：S 只能位于两端（与端点之间只允许 H），H 只能是最前/最后；
- 占位 CIGAR（kSmN，k==l_seq，m==65535）按 §4.2.2 用 CG:B,I 恢复；
- SEQ 保留 = 及歧义字母，奇数末半字节忽略；QUAL 仅收 0~93 或整列 255；
- 辅助字段支持 §4.2.4 全部类型（A,c,C,s,S,i,I,f,Z,H,B），同记录重标签拒绝。

出错抛 RecordError（.offset 为记录 vo）；头部错误抛 HeaderError（.offset 恒 0）。
"""

import bisect
import struct

__all__ = ["HeaderError", "RecordError", "parse_bam"]

# SEQ 半字节编码表（保留 = 与歧义字母）
_NT16 = "=ACMGRSVTWYHKDBN"

# CIGAR 操作码
_M, _I, _D, _N, _S, _H, _P, _EQ, _X = range(9)
_QUERY_OPS = (_M, _I, _S, _EQ, _X)    # 消费 SEQ
_REF_OPS = (_M, _D, _N, _EQ, _X)      # 消费参考
_ALIGN_OPS = (_M, _EQ, _X)            # 产生投影项

# B 型辅助数组的元素长度
_B_SIZES = {ord("c"): 1, ord("C"): 1, ord("s"): 2, ord("S"): 2,
            ord("i"): 4, ord("I"): 4, ord("f"): 4}

_TAG_CG = (ord("C") << 8) | ord("G")

_CORE = struct.Struct("<iiBBHHHiiii")  # 记录定长 32 字节核心
_I32 = struct.Struct("<i")


class HeaderError(Exception):
    """BAM 头部错误；按规定上报偏移 0。"""

    offset = 0

    def __init__(self):
        super().__init__("BAM header error")


class RecordError(Exception):
    """记录错误；offset 为出错记录的虚拟偏移 vo。"""

    def __init__(self, vo: int):
        super().__init__("BAM record error at virtual offset %d" % vo)
        self.offset = vo


def _is_alpha(b: int) -> bool:
    return 0x41 <= b <= 0x5A or 0x61 <= b <= 0x7A


def _is_alnum(b: int) -> bool:
    return _is_alpha(b) or 0x30 <= b <= 0x39


def _parse_header(data: bytes):
    """解析 BAM 二进制头，返回 (refs, pos)。refs 为 [(name, l_ref)]。"""
    n = len(data)
    if n < 8 or data[:4] != b"BAM\x01":
        raise HeaderError()
    (l_text,) = _I32.unpack_from(data, 4)
    if l_text < 0 or 8 + l_text > n:
        raise HeaderError()
    pos = 8 + l_text  # 文本头跳过
    if pos + 4 > n:
        raise HeaderError()
    (n_ref,) = _I32.unpack_from(data, pos)
    pos += 4
    if n_ref < 0:
        raise HeaderError()
    refs = []
    seen = set()
    for _ in range(n_ref):
        if pos + 4 > n:
            raise HeaderError()
        (l_name,) = _I32.unpack_from(data, pos)
        pos += 4
        if l_name < 1 or pos + l_name + 4 > n:
            raise HeaderError()
        name_raw = bytes(data[pos:pos + l_name])
        pos += l_name
        if name_raw[-1] != 0:  # 参考名须 NUL 结尾
            raise HeaderError()
        try:
            name = name_raw[:-1].decode("utf-8")
        except UnicodeDecodeError:
            raise HeaderError() from None
        if name in seen:  # 参考名唯一
            raise HeaderError()
        seen.add(name)
        (l_ref,) = _I32.unpack_from(data, pos)
        pos += 4
        if l_ref <= 0:  # 长度正
            raise HeaderError()
        refs.append((name, l_ref))
    return refs, pos


def _vo_locator(blocks):
    """由块表 [(ustart, usize, cstart)] 构造解压偏移 → vo 的换算函数。"""
    ustarts = [b[0] for b in blocks]

    def vo_of(u: int) -> int:
        i = bisect.bisect_right(ustarts, u) - 1
        ustart, _usize, cstart = blocks[i]
        return (cstart << 16) | (u - ustart)

    return vo_of


def _parse_aux(rec: bytes, p: int, vo: int):
    """解析辅助字段区 [p, len(rec))，校验类型与边界、拒绝重标签。

    返回 CG:B,I 的 uint32 数组（无则为 None）。
    """
    size = len(rec)
    tags = set()
    cg = None
    while p < size:
        if p + 3 > size:
            raise RecordError(vo)
        t0, t1, ty = rec[p], rec[p + 1], rec[p + 2]
        if not _is_alpha(t0) or not _is_alnum(t1):
            raise RecordError(vo)
        tag = (t0 << 8) | t1
        if tag in tags:  # 同记录重标签拒绝
            raise RecordError(vo)
        tags.add(tag)
        p += 3
        if ty == ord("A") or ty == ord("c") or ty == ord("C"):
            if p + 1 > size:
                raise RecordError(vo)
            p += 1
        elif ty == ord("s") or ty == ord("S"):
            if p + 2 > size:
                raise RecordError(vo)
            p += 2
        elif ty == ord("i") or ty == ord("I") or ty == ord("f"):
            if p + 4 > size:
                raise RecordError(vo)
            p += 4
        elif ty == ord("Z") or ty == ord("H"):  # NUL 结尾字符串
            q = rec.find(b"\0", p)
            if q < 0 or q >= size:
                raise RecordError(vo)
            p = q + 1
        elif ty == ord("B"):  # 数组：子类型 + int32 计数 + 元素
            if p + 5 > size:
                raise RecordError(vo)
            esz = _B_SIZES.get(rec[p])
            if esz is None:
                raise RecordError(vo)
            (count,) = _I32.unpack_from(rec, p + 1)
            if count < 0 or p + 5 + count * esz > size:
                raise RecordError(vo)
            if tag == _TAG_CG and rec[p] == ord("I"):
                cg = struct.unpack_from("<%dI" % count, rec, p + 5)
            p += 5 + count * esz
        else:
            raise RecordError(vo)
    return cg


def _check_clipping(ops, vo: int):
    """剪切位置遵循规范：H 只能是最前/最后操作，S 与端点之间只允许 H。"""
    n = len(ops)
    i = 0
    while i < n and ops[i][1] == _H:
        i += 1
    while i < n and ops[i][1] == _S:
        i += 1
    j = n
    while j > i and ops[j - 1][1] == _H:
        j -= 1
    while j > i and ops[j - 1][1] == _S:
        j -= 1
    for k in range(i, j):
        if ops[k][1] == _S or ops[k][1] == _H:
            raise RecordError(vo)


def _project(ops, pos0: int, l_ref: int, l_seq: int, seq_raw: bytes, vo: int):
    """校验 CIGAR 并把读段投影到参考位置，返回 columns = [[r, q, b], ...]。"""
    consume_q = 0
    span = 0
    decoded = []
    for w in ops:
        ln = w >> 4
        op = w & 0xF
        if ln == 0 or op > 8:  # 长度须正；B(9) 及更高操作码非法
            raise RecordError(vo)
        decoded.append((ln, op))
        if op in _QUERY_OPS:
            consume_q += ln
        if op in _REF_OPS:
            span += ln
    if consume_q != l_seq:  # 消费 SEQ 总长须一致
        raise RecordError(vo)
    if pos0 + span > l_ref:  # 参考跨度不得越界
        raise RecordError(vo)
    _check_clipping(decoded, vo)
    # 解码 SEQ（保留 = 与歧义字母；奇数末半字节忽略）
    seq = []
    for i in range(l_seq):
        byte = seq_raw[i >> 1]
        nib = (byte >> 4) if (i & 1) == 0 else (byte & 0xF)
        seq.append(_NT16[nib])
    # 投影：负链不再反补，q 为存储 SEQ 下标
    columns = []
    r = pos0
    q = 0
    for ln, op in decoded:
        if op in _ALIGN_OPS:  # M、=、X 产生项
            for _ in range(ln):
                columns.append([r, q, seq[q]])
                r += 1
                q += 1
        elif op == _I or op == _S:  # 仅推进 q
            q += ln
        elif op == _D or op == _N:  # 仅推进 r
            r += ln
        # H、P 均不推进
    return columns


def _parse_record(rec: bytes, refs, vo: int):
    """解析单条记录（block_size 之后的 block_size 字节）。

    跳过（返回 None）或产出 {"name", "ref", "vo", "columns"}。
    """
    (ref_id, pos0, l_read_name, _mapq, _bin, n_cigar, flag, l_seq,
     _next_ref, _next_pos, _tlen) = _CORE.unpack_from(rec, 0)
    if flag & 0x4:  # 未比对：跳过
        return None
    if l_seq == 0:  # SEQ 为空：跳过
        return None
    if l_seq < 0:
        raise RecordError(vo)
    if not 0 <= ref_id < len(refs):
        raise RecordError(vo)
    ref_name, l_ref = refs[ref_id]
    if pos0 < 0:
        raise RecordError(vo)
    size = len(rec)
    p = _CORE.size
    if l_read_name < 1 or p + l_read_name > size:
        raise RecordError(vo)
    name_raw = bytes(rec[p:p + l_read_name])
    p += l_read_name
    if name_raw[-1] != 0:  # 读名须 NUL 结尾
        raise RecordError(vo)
    try:
        name = name_raw[:-1].decode("utf-8")
    except UnicodeDecodeError:
        raise RecordError(vo) from None
    if p + 4 * n_cigar > size:
        raise RecordError(vo)
    ops = list(struct.unpack_from("<%dI" % n_cigar, rec, p)) if n_cigar else []
    p += 4 * n_cigar
    seq_nbytes = (l_seq + 1) // 2
    if p + seq_nbytes + l_seq > size:
        raise RecordError(vo)
    seq_raw = rec[p:p + seq_nbytes]
    p += seq_nbytes
    qual = rec[p:p + l_seq]
    p += l_seq
    cg = _parse_aux(rec, p, vo)
    # 占位 CIGAR（§4.2.2）：kSmN，k == l_seq，m == 65535 → 用 CG:B,I 恢复
    if (len(ops) == 2
            and ops[0] == ((l_seq << 4) | _S)
            and ops[1] == ((0xFFFF << 4) | _N)):
        if cg is None:
            raise RecordError(vo)
        ops = list(cg)
    # QUAL 仅收 0~93 或整列 255
    if not (all(b == 0xFF for b in qual) or all(b <= 93 for b in qual)):
        raise RecordError(vo)
    columns = _project(ops, pos0, l_ref, l_seq, seq_raw, vo)
    return {"name": name, "ref": ref_name, "vo": str(vo), "columns": columns}


def parse_bam(data: bytes, blocks):
    """解析解压后的完整 BAM 字节流，按文件序返回投影记录数组（可为 []）。"""
    refs, pos = _parse_header(data)
    vo_of = _vo_locator(blocks)
    out = []
    n = len(data)
    while pos < n:
        start = pos
        if n - start < 4:
            raise RecordError(vo_of(start))
        (block_size,) = _I32.unpack_from(data, start)
        if block_size < _CORE.size or start + 4 + block_size > n:
            raise RecordError(vo_of(start))
        rec = data[start + 4: start + 4 + block_size]
        vo = vo_of(start)
        item = _parse_record(rec, refs, vo)
        if item is not None:
            out.append(item)
        pos = start + 4 + block_size
    return out
