"""命令行与库入口：成功原子替换输出，失败保留旧文件并上报偏移。

用法：python -m bamproj <in.bam> <out.json>

退出码：
- 0：成功，输出已原子替换；
- 1：输入非法或 I/O 失败，stderr 打印十进制偏移
    （BGZF 块错 = 压缩文件零基块起点；记录错 = vo；头错 = 0），
    输出文件保持原样；
- 2：参数错误。
"""

import json
import os
import sys

from .bam import HeaderError, RecordError, parse_bam
from .bgzf import BgzfError, read_bgzf

__all__ = ["convert", "main"]


def convert(in_path: str, out_path: str) -> None:
    """把 BGZF 压缩 BAM 转为投影 JSON 并原子写入 out_path。

    失败抛 BgzfError / HeaderError / RecordError（均带 .offset），
    输出文件保持原样。
    """
    with open(in_path, "rb") as fh:
        raw = fh.read()
    data, blocks = read_bgzf(raw)
    records = parse_bam(data, blocks)
    payload = json.dumps(records, ensure_ascii=False, separators=(",", ":"))
    _atomic_write(out_path, payload.encode("utf-8"))


def _atomic_write(path: str, payload: bytes) -> None:
    """先写同目录临时文件并 fsync，再 os.replace 原子替换；失败清理临时文件。"""
    tmp = "%s.tmp.%d" % (path, os.getpid())
    try:
        with open(tmp, "wb") as fh:
            fh.write(payload)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, path)
        try:  # 尽量 fsync 目录，保证替换落盘
            dfd = os.open(os.path.dirname(os.path.abspath(path)) or ".", os.O_RDONLY)
            try:
                os.fsync(dfd)
            finally:
                os.close(dfd)
        except OSError:
            pass
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


def main(argv=None) -> int:
    argv = sys.argv if argv is None else argv
    if len(argv) != 3:
        print("usage: %s <in.bam> <out.json>" % os.path.basename(argv[0]),
              file=sys.stderr)
        return 2
    try:
        convert(argv[1], argv[2])
    except (BgzfError, HeaderError, RecordError) as exc:
        print(exc.offset, file=sys.stderr)
        return 1
    except OSError:
        print(0, file=sys.stderr)  # 输入不可读等 I/O 失败：按约定报 0
        return 1
    return 0
