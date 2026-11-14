"""旧参考 → 新参考的单样本 VCF 坐标转换（仅依赖标准库）。

处理流程：
1. 解析旧、新两份单条染色体参考 FASTA（纯 A/C/G/T，大小上限沿用核对工具）。
2. 解析制表符（TSV）坐标映射表：每行一个映射块，给出旧、新两侧的
   1 基闭区间与方向（+ / -）；两侧区间等长，同一侧各块互不重叠。
3. 解析旧参考上的单样本 VCF，沿用 REF 校验与未分相二倍体 GT 规则，
   仅接受每行单 ALT 的单碱基替换或事件长度 <= 20 bp 的纯插入 / 纯缺失。
4. 逐记录按映射块转换坐标：
   - 正向块：区间按位置一一对应；
   - 反向块：坐标镜像、等位反向互补，插缺在新链上重新锚定；
   转换后对照新参考做 REF 校验并左对齐归一化。
5. 同一目标归一变异合并剂量；合并后剂量超过二倍体上限 2 即整批失败。
6. 全部记录通过后才写出目标 VCF（带源行号 SRC 与 GT），不生成任何部分文件。

任何硬性错误都抛出 CheckError，代表整批转换失败。
"""

from dataclasses import dataclass

from .bio import (
    _BASES,
    CheckError,
    Reference,
    classify,
    normalize,
    parse_diploid_gt,
    parse_fasta,
)

# 映射表体积上限（字节），与单份 VCF 相同
MAX_MAP_BYTES = 256 * 1024

# 允许的映射表列名（可选表头，大小写不敏感）
_MAP_HEADER = ("old_start", "old_end", "new_start", "new_end", "strand")

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


# --------------------------------------------------------------------------- #
# 序列小工具
# --------------------------------------------------------------------------- #
def reverse_complement(seq: str) -> str:
    """纯 A/C/G/T 序列的反向互补。"""
    return seq.translate(_COMPLEMENT)[::-1]


# --------------------------------------------------------------------------- #
# 映射块
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Block:
    """一个坐标映射块；坐标均为 1 基闭区间。"""

    old_start: int
    old_end: int
    new_start: int
    new_end: int
    strand: str  # "+" | "-"
    line_no: int = 0

    @property
    def length(self) -> int:
        return self.old_end - self.old_start + 1

    def contains_old(self, start: int, end: int) -> bool:
        return self.old_start <= start and end <= self.old_end

    def map_forward(self, old_pos: int) -> int:
        """正向块内的单点坐标映射（1-based → 1-based）。"""
        return self.new_start + (old_pos - self.old_start)

    def map_reverse(self, old_pos: int) -> int:
        """反向块内的单点坐标镜像（1-based → 1-based）。"""
        return self.new_end - (old_pos - self.old_start)


def _parse_int_field(value: str, line_no: int, label: str) -> int:
    try:
        n = int(value)
    except ValueError:
        raise CheckError(
            f"映射表第 {line_no} 行：{label}={value!r} 不是整数。"
        )
    if n < 1:
        raise CheckError(f"映射表第 {line_no} 行：{label}={n} 不是正整数（1 基坐标）。")
    return n


def parse_mapping(data: bytes, old_ref: Reference, new_ref: Reference) -> list[Block]:
    """解析 TSV 坐标映射表并做完整性校验，返回按旧坐标升序排列的映射块。"""
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise CheckError(f"映射表不是 UTF-8 文本（字节 {exc.start} 处解码失败）。")

    blocks: list[Block] = []
    for line_no, raw in enumerate(text.splitlines(), start=1):
        line = raw.rstrip("\r")
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 5:
            raise CheckError(
                f"映射表第 {line_no} 行：应为 5 个制表符分隔列"
                f"（旧起点\\t旧终点\\t新起点\\t新终点\\t方向），实际 {len(fields)} 列。"
            )
        # 可选表头：首个数据行若列名完全匹配则跳过
        if not blocks and tuple(f.strip().lower() for f in fields) == _MAP_HEADER:
            continue

        old_start = _parse_int_field(fields[0].strip(), line_no, "旧起点")
        old_end = _parse_int_field(fields[1].strip(), line_no, "旧终点")
        new_start = _parse_int_field(fields[2].strip(), line_no, "新起点")
        new_end = _parse_int_field(fields[3].strip(), line_no, "新终点")
        strand = fields[4].strip()

        if strand not in ("+", "-"):
            raise CheckError(
                f"映射表第 {line_no} 行：方向 {strand!r} 非法，仅允许 '+'（正向）或 '-'（反向）。"
            )
        if old_end < old_start:
            raise CheckError(f"映射表第 {line_no} 行：旧区间终点 {old_end} 小于起点 {old_start}。")
        if new_end < new_start:
            raise CheckError(f"映射表第 {line_no} 行：新区间终点 {new_end} 小于起点 {new_start}。")

        old_len = old_end - old_start + 1
        new_len = new_end - new_start + 1
        if old_len != new_len:
            raise CheckError(
                f"映射表第 {line_no} 行：两侧区间等长要求不满足"
                f"（旧侧 {old_len} bp，新侧 {new_len} bp）。"
            )
        if old_end > len(old_ref.sequence):
            raise CheckError(
                f"映射表第 {line_no} 行：旧区间 [{old_start},{old_end}] 超出旧参考长度 "
                f"{len(old_ref.sequence)}。"
            )
        if new_end > len(new_ref.sequence):
            raise CheckError(
                f"映射表第 {line_no} 行：新区间 [{new_start},{new_end}] 超出新参考长度 "
                f"{len(new_ref.sequence)}。"
            )

        blocks.append(
            Block(old_start, old_end, new_start, new_end, strand, line_no=line_no)
        )

    if not blocks:
        raise CheckError("映射表不包含任何映射块。")

    blocks.sort(key=lambda b: b.old_start)

    # 旧侧互不重叠（端点相接允许）
    prev = blocks[0]
    for blk in blocks[1:]:
        if blk.old_start <= prev.old_end:
            raise CheckError(
                f"映射表第 {blk.line_no} 行与第 {prev.line_no} 行：旧侧区间重叠"
                f"（[{prev.old_start},{prev.old_end}] 与 [{blk.old_start},{blk.old_end}]）。"
            )
        prev = blk

    # 新侧块顺序可能因重排而与旧侧不同：按新区间起点排序后再查相邻重叠
    by_new = sorted(blocks, key=lambda b: b.new_start)
    prev = by_new[0]
    for blk in by_new[1:]:
        if blk.new_start <= prev.new_end:
            raise CheckError(
                f"映射表第 {blk.line_no} 行与第 {prev.line_no} 行：新侧区间重叠"
                f"（[{prev.new_start},{prev.new_end}] 与 "
                f"[{blk.new_start},{blk.new_end}]）。"
            )
        prev = blk

    return blocks


# --------------------------------------------------------------------------- #
# 源 VCF（旧参考，单样本、单 ALT）
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class SourceRecord:
    """源 VCF 中一条待转换记录（已经过 REF/GT 等基础校验）。"""

    line_no: int
    pos: int
    ref: str
    alt: str
    dosage: int
    gt: str


def parse_source_vcf(
    data: bytes, filename: str, refobj: Reference
) -> tuple[list[SourceRecord], str]:
    """解析旧参考上的单样本、单 ALT VCF，返回待转换记录与样本名。"""
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise CheckError(f"{filename} 不是 UTF-8 文本（字节 {exc.start} 处解码失败）。")

    lines = text.splitlines()
    header_idx = -1
    for i, line in enumerate(lines):
        if line.startswith("#CHROM"):
            header_idx = i
            break
    if header_idx < 0:
        raise CheckError(f"{filename}：缺少 '#CHROM' 表头行，不是合法 VCF。")

    header_fields = lines[header_idx].rstrip("\r").split("\t")
    if len(header_fields) < 10:
        raise CheckError(
            f"{filename}：'#CHROM' 表头只有 {len(header_fields)} 列，"
            "本工具要求恰好一个样本列（共 10 列）。"
        )
    if len(header_fields) > 10:
        raise CheckError(
            f"{filename}：'#CHROM' 表头有 {len(header_fields) - 9} 个样本列，"
            "本工具仅接受单样本 VCF。"
        )
    sample_name = header_fields[9]

    records: list[SourceRecord] = []
    for offset, raw in enumerate(lines[header_idx + 1:], start=header_idx + 2):
        line = raw.rstrip("\r")
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 8:
            raise CheckError(f"第 {offset} 行：列数不足 8 列，不是合法 VCF 记录。")

        chrom, pos_s, _vid, ref, alt_s, _qual, _filt, _info = fields[:8]
        if chrom != refobj.name:
            raise CheckError(
                f"第 {offset} 行：CHROM={chrom!r} 与旧参考 FASTA 名称 {refobj.name!r} 不一致。"
            )

        try:
            pos = int(pos_s)
        except ValueError:
            raise CheckError(f"第 {offset} 行：POS={pos_s!r} 不是整数。")
        if pos < 1 or pos > len(refobj.sequence):
            raise CheckError(
                f"第 {offset} 行：POS={pos} 超出旧参考范围（参考长度 {len(refobj.sequence)}）。"
            )

        if not ref or any(ch not in _BASES for ch in ref):
            raise CheckError(f"第 {offset} 行：REF 等位 {ref!r} 不是纯 A/C/G/T 序列。")
        alts = alt_s.split(",")
        if len(alts) != 1:
            raise CheckError(
                f"第 {offset} 行：本工具仅接受每行单 ALT，发现 {len(alts)} 个 ALT（{alt_s}）。"
            )
        alt = alts[0]
        if not alt or any(ch not in _BASES for ch in alt):
            raise CheckError(
                f"第 {offset} 行：ALT {alt_s!r} 为空或含非 A/C/G/T 字符"
                "（本工具不接受 *、<DEL>、符号化等位及缺失占位）。"
            )

        end = pos + len(ref) - 1
        if end > len(refobj.sequence):
            raise CheckError(
                f"第 {offset} 行：REF 延伸至 {end}，超出旧参考长度 {len(refobj.sequence)}。"
            )
        ref_sub = refobj.sequence[pos - 1:end]
        if ref_sub != ref:
            raise CheckError(
                f"第 {offset} 行：REF={ref} 与旧参考序列在 pos={pos} 处的 {ref_sub} 不一致，"
                "整批转换失败。"
            )

        if len(fields) == 8:
            raise CheckError(
                f"第 {offset} 行：缺少 FORMAT 与样本列，本工具要求单样本二倍体 GT。"
            )
        fmt_keys = fields[8].split(":")
        if not fmt_keys or fmt_keys[0] != "GT":
            raise CheckError(f"第 {offset} 行：FORMAT 必须以 GT 开头，实际为 {fields[8]!r}。")
        sample_values = fields[9].split(":")
        if len(sample_values) > len(fmt_keys):
            raise CheckError(f"第 {offset} 行：样本列的子字段多于 FORMAT 定义。")
        gt = sample_values[0] if sample_values else ""
        dosages = parse_diploid_gt(gt, offset, 1)
        dosage = dosages.get(1, 0)  # 单 ALT：要么剂量 0（0/0），要么 1/2

        records.append(SourceRecord(offset, pos, ref, alt, dosage, gt))

    return records, sample_name


# --------------------------------------------------------------------------- #
# 单记录坐标转换
# --------------------------------------------------------------------------- #
def _locate_block(blocks: list[Block], start: int, end: int, line_no: int) -> Block:
    """定位完整覆盖旧侧 [start, end] 的映射块；无映射或跨块即整批失败。"""
    for blk in blocks:
        if blk.old_start > end:
            break
        if blk.contains_old(start, end):
            return blk
    raise CheckError(
        f"第 {line_no} 行：变异旧参考区间 [{start},{end}] 未被任何映射块完整覆盖"
        "（无映射或变异跨块），整批转换失败。"
    )


def transform_record(
    record: SourceRecord,
    blocks: list[Block],
    new_ref: Reference,
):
    """转换一条源记录，返回归一后的 (Variant 三元组, kind, 剂量)。

    反向块：坐标镜像、等位反向互补，插缺在新链上重新锚定；
    随后统一对照新参考做 REF 校验与左对齐归一化。
    """
    start = record.pos
    end = record.pos + len(record.ref) - 1
    blk = _locate_block(blocks, start, end, record.line_no)
    seq = new_ref.sequence

    if blk.strand == "+":
        new_pos = blk.map_forward(record.pos)
        new_ref_allele = record.ref
        new_alt_allele = record.alt
    else:
        # 反向块：旧区间 [pos, end] 镜像到新链 [map(end), map(pos)]
        new_pos = blk.map_reverse(end)
        new_ref_allele = reverse_complement(record.ref)
        new_alt_allele = reverse_complement(record.alt)

    # 新参考边界与 REF 校验（先于归一化，给出目标 REF 不符的明确报错）
    new_end = new_pos + len(new_ref_allele) - 1
    if new_end > len(seq):
        raise CheckError(
            f"第 {record.line_no} 行：转换后 REF 延伸至 {new_end}，"
            f"超出新参考长度 {len(seq)}，整批转换失败。"
        )
    target_sub = seq[new_pos - 1:new_end]
    if target_sub != new_ref_allele:
        raise CheckError(
            f"第 {record.line_no} 行：目标 REF 不符——转换后在新参考 pos={new_pos} "
            f"处应为 {target_sub}，实际等位 REF={new_ref_allele}，整批转换失败。"
        )

    npos, nref, nalt = normalize(new_pos, new_ref_allele, new_alt_allele, seq)
    kind = classify(nref, nalt)  # 形态与 20 bp 长度上限校验
    return npos, nref, nalt, kind


# --------------------------------------------------------------------------- #
# 整批转换
# --------------------------------------------------------------------------- #
def _target_gt(dosage: int) -> str:
    if dosage == 1:
        return "0/1"
    return "1/1"


def convert_vcf(
    old_fasta: bytes,
    new_fasta: bytes,
    mapping_bytes: bytes,
    vcf_bytes: bytes,
    vcf_filename: str = "source.vcf",
) -> bytes:
    """整批转换入口；全部成功后才返回目标 VCF 字节，任何错误抛 CheckError。"""
    old_ref = parse_fasta(old_fasta)
    new_ref = parse_fasta(new_fasta)
    blocks = parse_mapping(mapping_bytes, old_ref, new_ref)
    records, sample_name = parse_source_vcf(vcf_bytes, vcf_filename, old_ref)

    # 目标变异 → (剂量, 源行号列表)
    merged: dict[tuple, list] = {}

    for record in records:
        if record.dosage == 0:
            # 0/0 记录不产生任何变异，与核对工具的剂量 0 忽略规则一致
            continue
        npos, nref, nalt, kind = transform_record(record, blocks, new_ref)
        key = (npos, nref, nalt, kind)
        entry = merged.get(key)
        if entry is None:
            merged[key] = [record.dosage, [record.line_no]]
            continue
        if entry[0] + record.dosage > 2:
            raise CheckError(
                f"第 {record.line_no} 行：与源第 {', '.join(map(str, entry[1]))} 行"
                f"转换为同一目标变异 {npos}:{nref}>{nalt}，合并剂量 "
                f"{entry[0]} + {record.dosage} 超过二倍体上限 2，整批转换失败。"
            )
        entry[0] += record.dosage
        entry[1].append(record.line_no)

    return _write_target_vcf(merged, new_ref, sample_name)


def _write_target_vcf(merged: dict[tuple, list], new_ref: Reference, sample_name: str) -> bytes:
    header = [
        "##fileformat=VCFv4.2",
        f"##reference={new_ref.name}",
        '##INFO=<ID=SRC,Number=.,Type=Integer,Description="源 VCF 物理行号（多个为合并来源）">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="未分相二倍体基因型">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_name,
    ]

    body: list[str] = []
    for (pos, ref, alt, kind), (dosage, line_numbers) in sorted(
        merged.items(), key=lambda kv: (kv[0][0], kv[0][1], kv[0][2])
    ):
        src = ",".join(str(n) for n in sorted(line_numbers))
        body.append(
            f"{new_ref.name}\t{pos}\t.\t{ref}\t{alt}\t.\t.\tSRC={src}\tGT\t{_target_gt(dosage)}"
        )

    return ("\n".join(header + body) + "\n").encode("utf-8")
