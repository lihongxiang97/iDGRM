from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from idgrm.io import read_expression, read_pairs, read_samples


class PairInputTests(unittest.TestCase):
    def test_reads_legacy_duplicate_headers(self):
        with TemporaryDirectory() as temp:
            path = Path(temp) / "pairs.tsv"
            path.write_text(
                "Duplicate 1\tLocation\tDuplicate 2\tLocation\n"
                "AT1G1\tChr1:1\tAT1G2\tChr1:2\n",
                encoding="utf-8",
            )
            pairs = read_pairs(path)
            self.assertEqual(pairs.loc[0, "gene1"], "AT1G1")
            self.assertEqual(pairs.loc[0, "gene2"], "AT1G2")

    def test_infers_parent_daughter_roles(self):
        with TemporaryDirectory() as temp:
            path = Path(temp) / "pairs.tsv"
            path.write_text(
                "Parent copy\tLocation\tDaughter copy\tLocation\n"
                "P\tChr1:1\tD\tChr1:2\n",
                encoding="utf-8",
            )
            pairs = read_pairs(path)
            self.assertEqual(pairs.loc[0, "role1"], "parent")
            self.assertEqual(pairs.loc[0, "role2"], "daughter")

    def test_rejects_missing_pair_gene_id(self):
        with TemporaryDirectory() as temp:
            path = Path(temp) / "pairs.tsv"
            path.write_text("gene1\tgene2\nG1\t\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "missing gene ID"):
                read_pairs(path)

    def test_rejects_nonnumeric_expression_value(self):
        with TemporaryDirectory() as temp:
            path = Path(temp) / "expression.tsv"
            path.write_text("gene_id\tleaf\nG1\tnot-a-number\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "nonnumeric value"):
                read_expression(path)

    def test_rejects_missing_sample_id(self):
        with TemporaryDirectory() as temp:
            path = Path(temp) / "samples.tsv"
            path.write_text("sample_id\ttissue\n\tleaf\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "missing sample ID"):
                read_samples(path, ["leaf_1"])

    def test_rejects_empty_pair_table(self):
        with TemporaryDirectory() as temp:
            path = Path(temp) / "pairs.tsv"
            path.write_text("gene1\tgene2\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "no duplicate-gene pairs"):
                read_pairs(path)

    def test_reads_headerless_legacy_pairs_without_losing_first_row(self):
        with TemporaryDirectory() as temp:
            path = Path(temp) / "Ath-b.segmental.pairs"
            path.write_text(
                "AT1G01030\tWGD\tAT2G46870\tbeta\n"
                "AT4G01500\tWGD\tAT2G46860\tbeta\n",
                encoding="utf-8",
            )
            pairs = read_pairs(path)
            self.assertEqual(len(pairs), 2)
            self.assertEqual(pairs.loc[0, "gene1"], "AT1G01030")
            self.assertEqual(pairs.loc[0, "gene2"], "AT2G46870")
            self.assertTrue(pairs.attrs["headerless"])

    def test_reads_ragged_parent_daughter_table_without_column_shift(self):
        with TemporaryDirectory() as temp:
            path = Path(temp) / "Ath.tandem-NG_parent-daughter.pairs"
            path.write_text(
                "Parent copy\tLocation\tDaughter copy\tLocation\tE-value\n"
                "AT1G02960\tAth-Chr1:666807\tAT1G02965\tAth-Chr1:671661\t8.5e-13\tparent_daughter\n",
                encoding="utf-8",
            )
            pairs = read_pairs(path)
            self.assertEqual(pairs.loc[0, "gene1"], "AT1G02960")
            self.assertEqual(pairs.loc[0, "gene2"], "AT1G02965")
            self.assertEqual(pairs.loc[0, "role1"], "parent")
            self.assertEqual(pairs.loc[0, "role2"], "daughter")

    def test_maps_transposed_parental_headers_to_daughter_parent_roles(self):
        with TemporaryDirectory() as temp:
            path = Path(temp) / "Ath.transposed.pairs"
            path.write_text(
                "Transposed\tLocation\tParental\tLocation\tE-value\n"
                "AT1G02030\tAth-Chr1:355124\tAT3G60580\tAth-Chr3:22393823\t6.9e-28\n",
                encoding="utf-8",
            )
            pairs = read_pairs(path)
            self.assertEqual(pairs.loc[0, "gene1"], "AT3G60580")
            self.assertEqual(pairs.loc[0, "gene2"], "AT1G02030")
            self.assertEqual(pairs.loc[0, "role1"], "parent")
            self.assertEqual(pairs.loc[0, "role2"], "daughter")

    def test_deduplicates_exact_pair_rows_and_records_warning(self):
        with TemporaryDirectory() as temp:
            path = Path(temp) / "pairs.tsv"
            path.write_text(
                "gene1\tgene2\nG1\tG2\nG1\tG2\n",
                encoding="utf-8",
            )
            pairs = read_pairs(path)
            self.assertEqual(len(pairs), 1)
            self.assertIn("Removed 1 exact duplicate", pairs.attrs["warnings"][0])


if __name__ == "__main__":
    unittest.main()
