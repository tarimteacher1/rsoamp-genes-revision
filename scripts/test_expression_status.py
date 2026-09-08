import csv
import tempfile
import unittest
from pathlib import Path

from finalize_expression_catalogue import matrix_subset, test_status


class ExpressionStatusTests(unittest.TestCase):
    def audit(self, retained=True, reason="RETAINED_TOTAL_COUNT_GE_10"):
        return {"gene_id": "gene1", "retained_for_DE": str(retained), "prefilter_status": reason}

    def test_retained(self):
        self.assertEqual(test_status(self.audit(), {"pvalue":"0.01", "padj":"0.04"}), "TESTED_WITH_ADJUSTED_P")

    def test_zero_and_low_count_are_not_model_missing(self):
        for reason, expected in [("FILTERED_ALL_ZERO_COUNTS", "NOT_TESTED_ALL_ZERO_COUNTS"),
                                 ("FILTERED_TOTAL_COUNT_LT_10", "NOT_TESTED_LOW_TOTAL_COUNT")]:
            self.assertEqual(test_status(self.audit(False, reason), None), expected)

    def test_independent_filtering_is_distinct(self):
        self.assertEqual(test_status(self.audit(), {"pvalue":"0.5", "padj":"NA"}), "INDEPENDENTLY_FILTERED_PADJ_UNAVAILABLE")

    def test_missing_p_not_assumed_low_count(self):
        self.assertEqual(test_status(self.audit(), {"pvalue":"NA", "padj":"NA"}), "P_VALUE_UNAVAILABLE")

    def test_unexplained_missing_or_extra_de_row_fails(self):
        with self.assertRaises(ValueError):
            test_status(self.audit(), None)
        with self.assertRaises(ValueError):
            test_status(self.audit(False, "FILTERED_ALL_ZERO_COUNTS"), {"pvalue":"0.5", "padj":"0.9"})

    def test_matrix_missing_is_not_zero_or_unquantified(self):
        member = {"member_id":"A", "gene_id":"gene1", "family":"Defensin", "final_subclass":"primary"}
        audit = {"gene1":self.audit(False, "FILTERED_ALL_ZERO_COUNTS")}
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"input.csv"
            path.write_text("gene_id,sample1\ngene_other,4\n")
            output = Path(directory)/"output.tsv"
            matrix_subset(path, [member], output, audit, True)
            with output.open() as handle:
                row = next(csv.DictReader(handle,delimiter="\t"))
            self.assertEqual(row["quantification_status"], "QUANTIFIED")
            self.assertEqual(row["matrix_status"], "NOT_AVAILABLE_PREFILTERED")
            self.assertEqual(row["sample1"], "")
            with self.assertRaises(ValueError):
                matrix_subset(path, [member], output, audit, False)


if __name__ == "__main__":
    unittest.main()
