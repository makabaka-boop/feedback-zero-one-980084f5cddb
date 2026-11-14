"""bamproj —— BGZF 压缩 BAM → 参考投影 JSON 后端。

自写 BGZF 分帧、BAM 解析与投影（SAM/BAM 1.6 第 4 节），
仅用标准库 zlib 做 DEFLATE 解压，不使用任何 BAM 读取器。

对外入口：
- convert(in_path, out_path)：转换并原子替换输出；失败抛错且保留旧文件。
- CLI：python -m bamproj <in.bam> <out.json>

错误约定（异常均带 .offset 属性）：
- bgzf.BgzfError：压缩文件零基块起点
- bam.RecordError：出错记录的虚拟偏移 vo
- bam.HeaderError：恒为 0
"""

from .bgzf import BgzfError, read_bgzf
from .bam import HeaderError, RecordError, parse_bam
from .cli import convert, main

__all__ = [
    "BgzfError",
    "HeaderError",
    "RecordError",
    "read_bgzf",
    "parse_bam",
    "convert",
    "main",
]

__version__ = "0.1.0"
