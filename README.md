# bamproj —— BGZF 压缩 BAM → 参考投影 JSON 后端

从零自写的测序入湖投影后端：读取 BGZF 压缩的 BAM，把每条已比对读段投影到
参考位置，按文件序输出 `{name, ref, vo, columns}` 数组。BGZF 分帧、BAM 解析、
投影全部自写（依 SAM/BAM 1.6 第 4 节），仅用标准库 `zlib` 做 DEFLATE 解压，
不使用任何 BAM 读取器。

## 用法

```bash
python -m bamproj <in.bam> <out.json>
# 或安装后：bamproj <in.bam> <out.json>
```

- 成功：退出码 0，输出**原子替换**（同目录临时文件 + fsync + `os.replace`）。
- 失败：退出码 1，stderr 打印一个十进制偏移，输出文件**保持原样**：
  - BGZF 块错 → 压缩文件零基块起点；
  - 记录错 → 该记录的虚拟偏移 vo；
  - 头错 → 0；
  - 多错任选其一。
- 参数错误：退出码 2。

库 API：`bamproj.convert(in_path, out_path)`，失败抛
`BgzfError / HeaderError / RecordError`（均带 `.offset`）。

## 输出格式

```json
[{"name":"r001","ref":"chr1","vo":"98","columns":[[0,0,"A"],[1,2,"G"]]}]
```

- 文件序数组；无记录为 `[]`；重名保留。
- `ref` 为参考名；`vo` 为十进制字符串：block_size 首字节所在块的压缩起点
  左移 16 位，或块内解压偏移。
- `columns` 每项 `[r, q, b]`：仅 M、=、X 产生项；`r` 零基参考位置，`q` 零基
  存储 SEQ 下标，`b` 为该字母；负链不再反补。

## 解析与校验规则

- **BGZF**：仅收 FLG=4、恰一个 BC 子字段（SLEN=2）；校验 BSIZE、CRC32、
  ISIZE（≤65536）；末块须空，内部空块继续读取；记录可跨块。
- **头**：文本头跳过；参考名唯一、长度正。
- **记录**：仅收已比对（FLAG 0x4 除外）且 SEQ 非空（l_seq>0）记录。
  - CIGAR 长度须正，操作码仅 0..8（MIDNSHP=X）；I、S 仅推进 q，D、N 仅
    推进 r，H、P 均不推进；消费 SEQ 总长须等于 l_seq；参考跨度不得越界；
    剪切位置遵循规范（S 只在两端、与端点间只允许 H，H 只在最前/最后）。
  - 占位 CIGAR（kSmN，k=l_seq，m=65535）按 §4.2.2 用 CG:B,I 恢复。
  - SEQ 保留 `=` 及歧义字母，奇数末半字节忽略；QUAL 仅收 0~93 或整列 255。
  - 辅助字段支持 §4.2.4 全部类型（A,c,C,s,S,i,I,f,Z,H,B）；同记录重标签拒绝。
- 非法输入**整次拒绝**。

## 目录结构

```
bamproj/            后端包
  bgzf.py           BGZF 分帧（自写，zlib 解压）
  bam.py            BAM 解析 + 参考投影
  cli.py            入口：原子替换、错误偏移上报
tests/
  make_fixtures.py  夹具生成器（独立的 BGZF/BAM 写出端）
  test_valid.py     合法投影对照 pysam；vo 独立走查
  test_invalid.py   非法输入拒绝与偏移验收
  test_cli.py       原子性、退出码、库 API
  fixtures/         生成的夹具与 manifest.json
```

## 验收

```bash
python3 tests/make_fixtures.py   # 生成夹具（幂等）
python3 -m pytest tests -q       # 运行全部验收
```

- 合法投影逐条对照 pysam（含经 `tell()` 取得的虚拟偏移）；
- vo 另由测试侧独立 BGZF/BAM 走查核对；
- 63 个夹具：合法（单块/跨块/内部空块/CG 占位/随机记录/空输出）与非法
  （BGZF 各字段、头部、记录各类错误）逐一验收偏移。

## 限制

- 整个解压流在内存中处理（适合 GB 级以下文件；更大文件可改为流式）。
- 参考序列视为线性；不校验文本头内容（按要求跳过）。
