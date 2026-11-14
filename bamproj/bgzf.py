"""BGZF 分帧读取（自写实现，仅用 zlib 做 DEFLATE 解压）。

依据 SAM/BAM 1.6 规范第 4 节的 BGZF 约定：
- 每块是合法 gzip 成员：ID1=0x1f、ID2=0x8b、CM=8，FLG 必须恰为 4（FEXTRA）；
- extra 字段中必须恰有一个 BC 子字段（SLEN=2），其数据为 BSIZE（块总字节数 - 1）；
- 解压后长度须等于 ISIZE（且 ≤ 65536），CRC32 须匹配；
- 末块必须为空块（ISIZE=0）；内部空块跳过继续读；
- BAM 记录可跨块（本模块把各块解压流拼接后统一交给上层）。

出错抛 BgzfError，offset 为出错块的压缩文件零基起点。
"""

import zlib

__all__ = ["BgzfError", "read_bgzf"]

# gzip 固定头 10 字节 + XLEN 2 字节 + 最小 extra（一个 BC 子字段 6 字节）
_MIN_HEADER = 18
# BGZF 块解压后最大长度
_MAX_ISIZE = 65536


class BgzfError(Exception):
    """BGZF 层错误；offset 为出错块的压缩文件零基起点。"""

    def __init__(self, offset: int):
        super().__init__("BGZF error at compressed offset %d" % offset)
        self.offset = offset


def _u16(buf: bytes, p: int) -> int:
    return buf[p] | (buf[p + 1] << 8)


def _u32(buf: bytes, p: int) -> int:
    return int.from_bytes(buf[p:p + 4], "little")


def _read_block(data: bytes, off: int):
    """解析 off 处的一个 BGZF 块，返回 (plain, total_size)。出错抛 BgzfError(off)。"""
    n = len(data)
    if n - off < _MIN_HEADER:
        raise BgzfError(off)
    if data[off] != 0x1F or data[off + 1] != 0x8B or data[off + 2] != 8:
        raise BgzfError(off)
    if data[off + 3] != 4:  # FLG 必须恰为 FEXTRA，其余标志位一律不收
        raise BgzfError(off)
    xlen = _u16(data, off + 10)
    hdr_end = off + 12 + xlen
    if hdr_end + 8 > n:
        raise BgzfError(off)
    # 解析 extra 子字段，必须恰有一个 BC（SLEN=2）
    bsize = None
    n_bc = 0
    p = off + 12
    while p < hdr_end:
        if p + 4 > hdr_end:
            raise BgzfError(off)
        si1, si2 = data[p], data[p + 1]
        slen = _u16(data, p + 2)
        if p + 4 + slen > hdr_end:
            raise BgzfError(off)
        if si1 == 0x42 and si2 == 0x43:  # 'B' 'C'
            n_bc += 1
            if slen != 2:
                raise BgzfError(off)
            bsize = _u16(data, p + 4)
        p += 4 + slen
    if n_bc != 1:
        raise BgzfError(off)
    total = bsize + 1
    if total < (hdr_end - off) + 8:  # BSIZE 连头部+尾部都装不下
        raise BgzfError(off)
    if off + total > n:  # 块被截断
        raise BgzfError(off)
    payload = data[hdr_end: off + total - 8]
    crc_want = _u32(data, off + total - 8)
    isize = _u32(data, off + total - 4)
    if isize > _MAX_ISIZE:
        raise BgzfError(off)
    decomp = zlib.decompressobj(-15)  # 裸 DEFLATE
    try:
        plain = decomp.decompress(payload) + decomp.flush()
    except zlib.error:
        raise BgzfError(off) from None
    if not decomp.eof or decomp.unused_data:  # 流须恰好终止于块尾
        raise BgzfError(off)
    if len(plain) != isize:
        raise BgzfError(off)
    if (zlib.crc32(plain) & 0xFFFFFFFF) != crc_want:
        raise BgzfError(off)
    return plain, total


def read_bgzf(data: bytes):
    """解析整个 BGZF 字节串。

    返回 (decompressed, blocks)：
    - decompressed：所有非空块解压流的拼接；
    - blocks：[(ustart, usize, cstart)]，仅含非空块，按解压流顺序排列，
      用于把解压偏移换算为虚拟偏移 vo = (cstart << 16) | (u - ustart)。
    """
    if len(data) == 0:
        raise BgzfError(0)
    blocks = []
    out = bytearray()
    off = 0
    cstart = 0
    last_empty = False
    while off < len(data):
        cstart = off
        plain, total = _read_block(data, off)
        if plain:
            blocks.append((len(out), len(plain), cstart))
            out += plain
            last_empty = False
        else:
            last_empty = True  # 内部空块：跳过继续读
        off += total
    if not last_empty:  # 末块须空
        raise BgzfError(cstart)
    return bytes(out), blocks
