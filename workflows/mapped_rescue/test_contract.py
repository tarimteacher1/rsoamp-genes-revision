import gzip
import tempfile
import unittest
from pathlib import Path

from run_stage import fragment_id, paired_records, orf_coordinates, overlaps, support


class Read:
    def __init__(self, name, read1=True, nh=1, blocks=None):
        self.query_name = name
        self.is_read1, self.is_read2 = read1, not read1
        self.is_reverse = read1
        self.is_unmapped = self.is_secondary = self.is_supplementary = self.is_qcfail = False
        self.nh, self.blocks = nh, blocks or [(10, 20)]

    def has_tag(self, tag):
        return self.nh is not None

    def get_tag(self, tag):
        return self.nh

    def get_blocks(self):
        return self.blocks


class Bam:
    references = ["Chr01"]

    def __init__(self, reads):
        self.reads = reads

    def fetch(self, *args):
        return iter(self.reads)


class Contract(unittest.TestCase):
    def pair(self, root, a, b):
        paths = [root / "a.gz", root / "b.gz"]
        for path, data in zip(paths, [a, b]):
            with gzip.open(path, "wb") as handle:
                handle.write(data)
        return paths

    def test_fragment_names(self):
        self.assertEqual(fragment_id(b"@SRR123.45/1 extra\n"), b"SRR123.45")
        self.assertEqual(fragment_id(b"@SRR123.45/2 extra\n"), b"SRR123.45")

    def test_good_pairs(self):
        with tempfile.TemporaryDirectory() as directory:
            paths = self.pair(Path(directory), b"@SRR1.1/1\nACTG\n+\nIIII\n", b"@SRR1.1/2\nCAGT\n+\nIIII\n")
            self.assertEqual(len(list(paired_records(*paths))), 1)

    def test_wrong_mates(self):
        with tempfile.TemporaryDirectory() as directory:
            paths = self.pair(Path(directory), b"@SRR1.1\nACTG\n+\nIIII\n", b"@SRR1.2\nCAGT\n+\nIIII\n")
            with self.assertRaises(ValueError):
                list(paired_records(*paths))

    def test_truncated_and_quality(self):
        with tempfile.TemporaryDirectory() as directory:
            for bad in [b"@SRR1.1\nACTG\n+\n", b"@SRR1.1\nACTG\n+\nIII\n", b""]:
                paths = self.pair(Path(directory), b"@SRR1.1\nACTG\n+\nIIII\n", bad)
                with self.assertRaises(ValueError):
                    list(paired_records(*paths))

    def test_orf_coordinates(self):
        self.assertEqual(orf_coordinates("M6MERGE.1.1_2 [4 - 96]"), (4, 96))
        with self.assertRaises(ValueError):
            orf_coordinates("x [96 - 4] (REVERSE SENSE)")

    def test_half_open_overlap(self):
        self.assertFalse(overlaps([(10, 20)], [(20, 30)]))
        self.assertTrue(overlaps([(10, 21)], [(20, 30)]))

    def test_mate_dedup_and_multi(self):
        bam = Bam([Read("pairA"), Read("pairA", read1=False), Read("pairB", nh=4)])
        value = support(bam, "Chr01", [(10, 20)], "+")
        self.assertEqual(value["overlapping_fragments"], 2)
        self.assertEqual(value["unique_fragments"], 1)
        self.assertEqual(value["multimapping_fragments"], 1)
        self.assertEqual(value["unique_exon_coverage_fraction"], 1)

    def test_missing_nh_stops(self):
        with self.assertRaises(ValueError):
            support(Bam([Read("pairA", nh=None)]), "Chr01", [(10, 20)], "+")

    def test_spliced_blocks_not_introns(self):
        value = support(Bam([Read("pairA", blocks=[(10, 20), (100, 110)])]), "Chr01", [(30, 40)], "+")
        self.assertEqual(value["overlapping_fragments"], 0)

    def test_wrong_strand_excluded(self):
        value = support(Bam([Read("pairA")]), "Chr01", [(10, 20)], "-")
        self.assertEqual(value["overlapping_fragments"], 0)


if __name__ == "__main__":
    unittest.main()
