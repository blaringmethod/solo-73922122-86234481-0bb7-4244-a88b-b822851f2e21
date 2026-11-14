"""参考坐标转换的核心逻辑测试，仅依赖标准库：
python -I -m unittest -v tests.test_convert
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from app.bio import CheckError  # noqa: E402
from app.convert import (  # noqa: E402
    MAX_MAP_BYTES,
    convert_vcf,
    parse_mapping,
    parse_source_vcf,
    reverse_complement,
    transform_record,
)
from app.bio import parse_fasta  # noqa: E402

# 旧参考 1..8：A A C G T A C G
OLD = b">chrOld\nAACGTACG\n"
# 块 A：旧 [1,4] AACG 反向 -> 新 [1,4] = RC(AACG) = CGTT
# 块 B：旧 [5,8] TACG 正向 -> 新 [5,8] = TACG
NEW = b">chrNew\nCGTTTACG\n"
MAP = b"1\t4\t1\t4\t-\n5\t8\t5\t8\t+\n"

HEADER = (
    "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
)


def vcf(rows):
    return (HEADER + "".join(r + "\n" for r in rows)).encode()


def rec(pos, ref, alt, gt):
    return f"chrOld\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT\t{gt}"


class Fixture:
    def __init__(self):
        self.old = parse_fasta(OLD)
        self.new = parse_fasta(NEW)
        self.blocks = parse_mapping(MAP, self.old, self.new)


class ReverseComplementTests(unittest.TestCase):
    def test_basic(self):
        self.assertEqual(reverse_complement("AACG"), "CGTT")
        self.assertEqual(reverse_complement("A"), "T")
        self.assertEqual(reverse_complement("ATCG"), "CGAT")


class MappingTests(unittest.TestCase):
    def setUp(self):
        self.fx = Fixture()

    def test_valid_blocks_sorted(self):
        blocks = parse_mapping(
            b"5\t8\t5\t8\t+\n1\t4\t1\t4\t-\n", self.fx.old, self.fx.new
        )
        self.assertEqual([b.old_start for b in blocks], [1, 5])
        self.assertEqual(blocks[0].strand, "-")
        self.assertEqual(blocks[0].length, 4)

    def test_header_row_skipped(self):
        data = (
            b"old_start\told_end\tnew_start\tnew_end\tstrand\n" + MAP
        )
        blocks = parse_mapping(data, self.fx.old, self.fx.new)
        self.assertEqual(len(blocks), 2)

    def test_comment_and_blank_lines_skipped(self):
        data = b"# comment\n\n" + MAP
        self.assertEqual(len(parse_mapping(data, self.fx.old, self.fx.new)), 2)

    def test_wrong_column_count(self):
        with self.assertRaises(CheckError) as ctx:
            parse_mapping(b"1\t4\t1\t4\n", self.fx.old, self.fx.new)
        self.assertIn("5 个", str(ctx.exception))

    def test_bad_strand(self):
        with self.assertRaises(CheckError):
            parse_mapping(b"1\t4\t1\t4\t*\n", self.fx.old, self.fx.new)

    def test_unequal_length(self):
        with self.assertRaises(CheckError) as ctx:
            parse_mapping(b"1\t4\t1\t5\t-\n", self.fx.old, self.fx.new)
        self.assertIn("等长", str(ctx.exception))

    def test_old_side_overlap(self):
        with self.assertRaises(CheckError) as ctx:
            parse_mapping(
                b"1\t4\t1\t4\t+\n3\t6\t5\t8\t+\n", self.fx.old, self.fx.new
            )
        self.assertIn("旧侧区间重叠", str(ctx.exception))

    def test_new_side_overlap_after_reorder(self):
        # 旧侧不重叠，但新侧重排后区间重叠
        with self.assertRaises(CheckError) as ctx:
            parse_mapping(
                b"1\t4\t1\t4\t+\n5\t8\t3\t6\t+\n", self.fx.old, self.fx.new
            )
        self.assertIn("新侧区间重叠", str(ctx.exception))

    def test_adjant_blocks_touching_allowed(self):
        blocks = parse_mapping(
            b"1\t4\t1\t4\t+\n5\t8\t5\t8\t+\n", self.fx.old, self.fx.new
        )
        self.assertEqual(len(blocks), 2)

    def test_non_integer_coordinate(self):
        with self.assertRaises(CheckError):
            parse_mapping(b"x\t4\t1\t4\t+\n", self.fx.old, self.fx.new)

    def test_non_positive_coordinate(self):
        with self.assertRaises(CheckError):
            parse_mapping(b"0\t4\t1\t4\t+\n", self.fx.old, self.fx.new)

    def test_old_range_exceeds_reference(self):
        with self.assertRaises(CheckError):
            parse_mapping(b"1\t9\t1\t9\t+\n", self.fx.old, self.fx.new)

    def test_new_range_exceeds_reference(self):
        with self.assertRaises(CheckError):
            parse_mapping(b"1\t4\t6\t9\t+\n", self.fx.old, self.fx.new)

    def test_empty_mapping(self):
        with self.assertRaises(CheckError):
            parse_mapping(b"# only comment\n", self.fx.old, self.fx.new)

    def test_reversed_interval_rejected(self):
        with self.assertRaises(CheckError):
            parse_mapping(b"4\t1\t1\t4\t+\n", self.fx.old, self.fx.new)


class TransformTests(unittest.TestCase):
    def setUp(self):
        self.fx = Fixture()

    def _one(self, row):
        recs, _ = parse_source_vcf(vcf([row]), "x.vcf", self.fx.old)
        return transform_record(recs[0], self.fx.blocks, self.fx.new)

    def test_reverse_snv_mirrored_and_complemented(self):
        # 旧 pos3 C>G -> 镜像新 pos2，G 的互补为 C
        pos, ref, alt, kind = self._one(rec(3, "C", "G", "0/1"))
        self.assertEqual((pos, ref, alt), (2, "G", "C"))
        self.assertEqual(kind, "SNV")

    def test_forward_snv_kept(self):
        pos, ref, alt, kind = self._one(rec(5, "T", "A", "1/1"))
        self.assertEqual((pos, ref, alt), (5, "T", "A"))
        self.assertEqual(kind, "SNV")

    def test_reverse_insertion_reanchored(self):
        # 旧 pos4 G>GC（插 C）-> RC 后新 pos1 C>GC（插 G）
        pos, ref, alt, kind = self._one(rec(4, "G", "GC", "0/1"))
        self.assertEqual((pos, ref, alt), (1, "C", "GC"))
        self.assertEqual(kind, "INS")

    def test_forward_insertion_kept(self):
        # 旧 pos5 T>TA -> 新 pos5 T>TA
        pos, ref, alt, kind = self._one(rec(5, "T", "TA", "0/1"))
        self.assertEqual((pos, ref, alt), (5, "T", "TA"))
        self.assertEqual(kind, "INS")

    def test_reverse_deletion_reanchored(self):
        # 旧 pos3 CG>C（删 old4 的 G）-> 新 pos1 CG>G
        pos, ref, alt, kind = self._one(rec(3, "CG", "C", "0/1"))
        self.assertEqual((pos, ref, alt), (1, "CG", "G"))
        self.assertEqual(kind, "DEL")

    def test_forward_deletion_kept(self):
        pos, ref, alt, kind = self._one(rec(5, "TA", "T", "0/1"))
        self.assertEqual((pos, ref, alt), (5, "TA", "T"))
        self.assertEqual(kind, "DEL")

    def test_variant_spanning_blocks_fails(self):
        # 缺失区间 [4,5] 横跨两个块
        with self.assertRaises(CheckError) as ctx:
            self._one(rec(4, "GT", "G", "0/1"))
        self.assertIn("跨块", str(ctx.exception))

    def test_unmapped_position_fails(self):
        blocks = parse_mapping(
            b"1\t3\t2\t4\t-\n5\t8\t5\t8\t+\n", self.fx.old, self.fx.new
        )
        recs, _ = parse_source_vcf(vcf([rec(4, "G", "C", "0/1")]), "x.vcf", self.fx.old)
        with self.assertRaises(CheckError) as ctx:
            transform_record(recs[0], blocks, self.fx.new)
        self.assertIn("无映射", str(ctx.exception))

    def test_target_ref_mismatch_fails(self):
        # 新参考 pos2 由 G 改为 A：反向 SNV 转换后 REF=G 与新参考不符
        bad_new = parse_fasta(b">chrNew\nCATT TACG\n".replace(b" ", b""))
        with self.assertRaises(CheckError) as ctx:
            self._one_via(rec(3, "C", "G", "0/1"), bad_new)
        self.assertIn("目标 REF 不符", str(ctx.exception))

    def _one_via(self, row, new_ref):
        recs, _ = parse_source_vcf(vcf([row]), "x.vcf", self.fx.old)
        return transform_record(recs[0], self.fx.blocks, new_ref)

    def test_insertion_normalized_against_new_reference(self):
        # 新参考 pos3-5 为 TTT 串（CGTTTACG）。旧 pos1 A>AA（插 A）经反向块：
        # 镜像到新 pos4，RC 得 T>TT（插 T）；左对齐应越过 TTT 串，
        # 锚定到串前的 G（pos2），归一为 2:G>GT（与 bcftools/vt 左对齐一致）。
        pos, ref, alt, kind = self._one(rec(1, "A", "AA", "0/1"))
        self.assertEqual(kind, "INS")
        self.assertEqual((pos, ref, alt), (2, "G", "GT"))

    def test_reverse_deletion_left_aligned_inside_new_repeat(self):
        # 反向块旧 [1,4]=AACG -> 新 [1,4]=CGTT；新参考 pos3-5 为 TTT 串。
        # 旧 pos4 G>GT（在旧 ref pos4 后插 T），RC 后为新 pos1 C>CA？
        # 改为构造缺失：旧 pos3 CG>C（删 old4=G）已在他例验证；
        # 这里用旧 pos1 AA>A（删 old1=A），RC 后为新 pos4 TT>T（删一个 T），
        # 落在新参考 TTT 串上，归一化左对齐到串起点锚点 pos2：2:GT>G。
        pos, ref, alt, kind = self._one(rec(1, "AA", "A", "0/1"))
        self.assertEqual(kind, "DEL")
        self.assertEqual((pos, ref), (2, "GT"))
        self.assertEqual(alt, "G")

    def test_snv_not_moved_by_normalization(self):
        # SNV 形态归一化后位置与等位不变（反向链已完成镜像与互补）
        pos, ref, alt, kind = self._one(rec(6, "A", "C", "0/1"))
        self.assertEqual(kind, "SNV")
        self.assertEqual((pos, ref, alt), (6, "A", "C"))

    def test_indel_too_long_rejected(self):
        long_ins = "G" + "C" * 21
        with self.assertRaises(CheckError):
            self._one(rec(5, "T", "T" + "A" * 21, "0/1"))
        self.assertEqual(len(long_ins), 22)  # 事件 21 bp


class SourceVcfTests(unittest.TestCase):
    def setUp(self):
        self.fx = Fixture()

    def test_single_alt_required(self):
        with self.assertRaises(CheckError) as ctx:
            parse_source_vcf(vcf([rec(3, "C", "G,T", "1/2")]), "x.vcf", self.fx.old)
        self.assertIn("单 ALT", str(ctx.exception))

    def test_phased_gt_rejected(self):
        with self.assertRaises(CheckError):
            parse_source_vcf(vcf([rec(3, "C", "G", "0|1")]), "x.vcf", self.fx.old)

    def test_missing_gt_rejected(self):
        with self.assertRaises(CheckError):
            parse_source_vcf(vcf([rec(3, "C", "G", "./1")]), "x.vcf", self.fx.old)

    def test_ref_check_against_old_reference(self):
        with self.assertRaises(CheckError) as ctx:
            parse_source_vcf(vcf([rec(3, "G", "T", "0/1")]), "x.vcf", self.fx.old)
        self.assertIn("REF", str(ctx.exception))

    def test_zero_dose_record_kept_but_marked(self):
        recs, name = parse_source_vcf(vcf([rec(3, "C", "G", "0/0")]), "x.vcf", self.fx.old)
        self.assertEqual(name, "S1")
        self.assertEqual(recs[0].dosage, 0)

    def test_homozygous_dosage_two(self):
        recs, _ = parse_source_vcf(vcf([rec(3, "C", "G", "1/1")]), "x.vcf", self.fx.old)
        self.assertEqual(recs[0].dosage, 2)
        self.assertEqual(recs[0].gt, "1/1")

    def test_physical_line_number_retained(self):
        data = vcf([rec(3, "C", "G", "0/1")])
        recs, _ = parse_source_vcf(data, "x.vcf", self.fx.old)
        self.assertEqual(recs[0].line_no, 3)

    def test_multi_sample_rejected(self):
        data = (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
            + rec(3, "C", "G", "0/1") + "\n"
        ).encode()
        with self.assertRaises(CheckError):
            parse_source_vcf(data, "x.vcf", self.fx.old)


class ConvertBatchTests(unittest.TestCase):
    def test_happy_path_target_vcf_content(self):
        data = vcf([
            rec(3, "C", "G", "0/1"),     # -> 2:G>C dose1, SRC=3
            rec(5, "T", "A", "1/1"),     # -> 5:T>A dose2, SRC=4
            rec(2, "A", "G", "0/0"),     # dose0 忽略
        ])
        out = convert_vcf(OLD, NEW, MAP, data).decode()
        self.assertIn("##reference=chrNew", out)
        self.assertIn("##INFO=<ID=SRC", out)
        body = [ln for ln in out.splitlines() if ln and not ln.startswith("#")]
        self.assertEqual(len(body), 2)
        self.assertEqual(body[0], "chrNew\t2\t.\tG\tC\t.\t.\tSRC=3\tGT\t0/1")
        self.assertEqual(body[1], "chrNew\t5\t.\tT\tA\t.\t.\tSRC=4\tGT\t1/1")

    def test_records_sorted_by_position(self):
        data = vcf([
            rec(5, "T", "A", "0/1"),
            rec(3, "C", "G", "0/1"),
        ])
        out = convert_vcf(OLD, NEW, MAP, data).decode()
        body = [ln for ln in out.splitlines() if ln and not ln.startswith("#")]
        self.assertTrue(body[0].startswith("chrNew\t2\t"))
        self.assertTrue(body[1].startswith("chrNew\t5\t"))

    def test_merge_same_target_variant_dosage_and_src(self):
        # 两条源记录转换为同一目标变异：剂量 1+1=2，SRC 合并两个行号
        data = vcf([
            rec(3, "C", "G", "0/1"),
            rec(3, "C", "G", "0/1"),
        ])
        out = convert_vcf(OLD, NEW, MAP, data).decode()
        body = [ln for ln in out.splitlines() if ln and not ln.startswith("#")]
        self.assertEqual(len(body), 1)
        self.assertEqual(body[0], "chrNew\t2\t.\tG\tC\t.\t.\tSRC=3,4\tGT\t1/1")

    def test_merged_dosage_over_two_fails_entire_batch(self):
        data = vcf([
            rec(3, "C", "G", "1/1"),
            rec(3, "C", "G", "0/1"),
        ])
        with self.assertRaises(CheckError) as ctx:
            convert_vcf(OLD, NEW, MAP, data)
        self.assertIn("超过二倍体上限 2", str(ctx.exception))

    def test_cross_block_fails_entire_batch(self):
        data = vcf([rec(4, "GT", "G", "0/1")])
        with self.assertRaises(CheckError):
            convert_vcf(OLD, NEW, MAP, data)

    def test_target_ref_mismatch_fails_entire_batch(self):
        bad_new = b">chrNew\nCATTTACG\n"  # pos2 G->A
        data = vcf([rec(3, "C", "G", "0/1")])
        with self.assertRaises(CheckError):
            convert_vcf(OLD, bad_new, MAP, data)

    def test_empty_variant_set_emits_header_only(self):
        out = convert_vcf(OLD, NEW, MAP, vcf([rec(3, "C", "G", "0/0")])).decode()
        body = [ln for ln in out.splitlines() if ln and not ln.startswith("#")]
        self.assertEqual(body, [])
        self.assertIn("#CHROM", out)

    def test_map_size_constant_matches_vcf_limit(self):
        self.assertEqual(MAX_MAP_BYTES, 256 * 1024)


if __name__ == "__main__":
    unittest.main()
