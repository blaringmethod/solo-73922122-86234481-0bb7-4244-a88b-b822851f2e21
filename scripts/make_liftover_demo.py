"""生成参考坐标转换的演示数据（scripts/make_demo.py 的姊妹脚本）。

运行：python3 scripts/make_liftover_demo.py
产物（均写入 samples/）：
  old_ref.fa        旧参考（单染色体）
  new_ref.fa        新参考（单染色体；前段反向互补、后段正向重排）
  liftover.map.tsv  5 列 TSV 映射块：旧起点 旧终点 新起点 新终点 方向
  source.vcf        按旧参考调用的单样本、单 ALT VCF（覆盖 SNV/INS/DEL 与正反链）
所有源 VCF 记录的 REF 均对照旧参考自动计算，保证合法。
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from app.bio import parse_fasta  # noqa: E402
from app.convert import reverse_complement  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "samples"

# 旧参考拆成两个片段：A 段反向映射，B 段正向映射；新参考中 B 段排在 A 段之前，
# 以此演示重排 + 反向互补。片段内保留重复区，便于演示归一化。
SEG_A = "AAAACGTACGTACGTT"          # 16 bp，含 A 串
SEG_B = "GGATCCATATATATGCGC"        # 18 bp，含 AT 二碱基重复区
LEN_A, LEN_B = len(SEG_A), len(SEG_B)
OLD_SEQ = SEG_A + SEG_B             # 34 bp
# 反向块要求新链序列 = 旧 A 段的反向互补；正向块 B 段序列保持不变
NEW_SEQ = SEG_B + reverse_complement(SEG_A)

HEADER = (
    "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tDEMO_SAMPLE\n"
)


def row(pos, ref, alt, gt):
    return f"chrOld\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}"


def main() -> None:
    # 自洽性检查：按映射规则拼出的新参考必须与文件一致
    assert NEW_SEQ == SEG_B + reverse_complement(SEG_A)
    parse_fasta(">chrOld\n".encode() + OLD_SEQ.encode())  # 纯 ACGT 校验

    (OUT / "old_ref.fa").write_text(f">chrOld 转换演示旧参考\n{OLD_SEQ}\n", encoding="utf-8")
    (OUT / "new_ref.fa").write_text(f">chrNew 转换演示新参考（重排+反向互补）\n{NEW_SEQ}\n", encoding="utf-8")
    (OUT / "liftover.map.tsv").write_text(
        "old_start\told_end\tnew_start\tnew_end\tstrand\n"
        f"1\t{LEN_A}\t{LEN_B + 1}\t{LEN_B + LEN_A}\t-\n"
        f"{LEN_A + 1}\t{LEN_A + LEN_B}\t1\t{LEN_B}\t+\n",
        encoding="utf-8",
    )

    records = [
        # 反向块（A 段，旧 1..16 -> 新 19..34，方向 -）
        row(2, SEG_A[1], "G", "0/1"),                       # 反向 SNV
        row(4, SEG_A[3], SEG_A[3] + "T", "0/1"),            # 反向插入（插 T）
        row(8, SEG_A[7:9], SEG_A[7], "0/1"),                # 反向缺失（删 1 bp）
        # 正向块（B 段，旧 17..34 -> 新 1..18，方向 +）
        row(LEN_A + 2, SEG_B[1], "A", "1/1"),               # 正向 SNV，纯合
        row(LEN_A + 4, SEG_B[3], SEG_B[3] + "G", "0/1"),    # 正向插入（插 G）
        row(LEN_A + 9, SEG_B[8:10], SEG_B[8], "0/1"),       # 正向缺失（AT 重复区内）
    ]
    (OUT / "source.vcf").write_text(HEADER + "\n".join(records) + "\n", encoding="utf-8")
    print("wrote samples/old_ref.fa, samples/new_ref.fa, samples/liftover.map.tsv, samples/source.vcf")


if __name__ == "__main__":
    main()
