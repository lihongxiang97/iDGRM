from pathlib import Path
import tempfile
import unittest

import pandas as pd
from idgrm.counts import validate_count_inputs,fit_copy_counts
from idgrm.cli import build_parser


class CountInputsTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
        self.counts=self.root/"counts.tsv";self.pairs=self.root/"pairs.tsv";self.samples=self.root/"samples.tsv"
        pd.DataFrame([[10,12],[5,4],[30,31]],index=["a","b","bg"],columns=["s1","s2"]).to_csv(self.counts,sep="\t")
        pd.DataFrame(dict(pair_id=["p1"],gene1=["a"],gene2=["b"])).to_csv(self.pairs,sep="\t",index=False)
        pd.DataFrame(dict(sample_id=["s1","s2"],tissue=["leaf","leaf"],replicate=["1","2"])).to_csv(self.samples,sep="\t",index=False)
    def tearDown(self): self.tmp.cleanup()
    def check(self):return validate_count_inputs(self.counts,self.pairs,self.samples)
    def test_valid(self):self.assertEqual(self.check()["libraries"],2)
    def test_fractional_rejected(self):
        frame=pd.read_csv(self.counts,sep="\t",index_col=0);frame=frame.astype(float);frame.iloc[0,0]=1.2;frame.to_csv(self.counts,sep="\t")
        with self.assertRaisesRegex(ValueError,"integer"):self.check()
    def test_negative_rejected(self):
        frame=pd.read_csv(self.counts,sep="\t",index_col=0);frame.iloc[0,0]=-1;frame.to_csv(self.counts,sep="\t")
        with self.assertRaises(ValueError):self.check()
    def test_order_rejected(self):
        frame=pd.read_csv(self.samples,sep="\t");frame.iloc[::-1].to_csv(self.samples,sep="\t",index=False)
        with self.assertRaisesRegex(ValueError,"ordered"):self.check()
    def test_reciprocal_requires_two_matrices(self):
        with self.assertRaisesRegex(ValueError,"requires"):
            fit_copy_counts(self.counts,self.pairs,self.samples,self.root/"out",observability_method="reciprocal-common-cds")
    def test_path_like_tissue_rejected(self):
        frame=pd.read_csv(self.samples,sep="\t");frame["tissue"]="../bad";frame.to_csv(self.samples,sep="\t",index=False)
        with self.assertRaisesRegex(ValueError,"Tissue"):self.check()
    def test_replicate_reuse_rejected(self):
        frame=pd.read_csv(self.samples,sep="\t");frame["replicate"]=1;frame.to_csv(self.samples,sep="\t",index=False)
        with self.assertRaisesRegex(ValueError,"Replicate"):self.check()
    def test_counts_cli_separate_engine(self):
        args=build_parser().parse_args(["fit-copy-counts","--counts","c","--duplicate-pairs","p","--sample-metadata","s",
            "--observability-method","uncorrected","--output-dir","o"])
        self.assertEqual(args.test,"wald")


if __name__=="__main__":unittest.main()
