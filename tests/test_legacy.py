from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from idgrm.config import IDGRMConfig
from idgrm.io import read_pairs
from idgrm.legacy import load_legacy_deseq2
from idgrm.pipeline import run_legacy_analysis


class LegacyTests(unittest.TestCase):
    def test_imports_deprecated_deseq2_shape_without_splitting_gene_ids(self):
        with TemporaryDirectory() as temp_name:
            temp = Path(temp_name)
            pair_path = temp / "pairs.tsv"
            pair_path.write_text("gene1\tgene2\nG-1\tG-2\n", encoding="utf-8")
            header = '"",baseMean,log2FoldChange,lfcSE,stat,pvalue,padj\n'
            (temp / "Ath.td.pairs.root.DESeq2.csv").write_text(
                header + '"G-1-G-2",10,-2,0.1,-20,1e-8,1e-7\n', encoding="utf-8"
            )
            (temp / "Ath.td.pairs.leaf.DESeq2.csv").write_text(
                header + '"G-1-G-2",10,2,0.1,20,1e-8,1e-7\n', encoding="utf-8"
            )
            pairs = read_pairs(pair_path)
            evidence, warnings = load_legacy_deseq2(temp, pairs, IDGRMConfig())
            self.assertFalse(warnings)
            self.assertEqual(set(evidence["direction"]), {"gene1_high", "gene2_high"})
            self.assertEqual(set(evidence["pair_id"]), {"G-1__G-2"})

            result = run_legacy_analysis(
                deseq_directory=temp,
                pairs_path=pair_path,
                output_dir=temp / "out",
                config=IDGRMConfig(),
            )
            self.assertEqual(result.classifications.iloc[0]["science_class"], "SUB_OR_NEO")
            self.assertEqual(
                result.classifications.iloc[0]["extended_fate"], "AMBIGUOUS_SUB_NEO"
            )
            for output_path in result.output_paths.values():
                self.assertTrue(output_path.exists(), output_path)


if __name__ == "__main__":
    unittest.main()
