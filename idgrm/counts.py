"""Raw integer count fitting, distinct from the log-expression matrix path."""
from __future__ import annotations

import hashlib
import csv
import json
import os
from pathlib import Path
import re
import subprocess

import numpy as np
import pandas as pd


def validate_count_inputs(counts_path,pairs_path,samples_path,copy1_path=None,copy2_path=None):
    for path in (counts_path,copy1_path,copy2_path):
        if path:
            with Path(path).open(newline="",encoding="utf-8-sig") as handle:
                header=next(csv.reader(handle,delimiter="\t"))
            if len(header)!=len(set(header)):raise ValueError("Duplicate matrix column labels")
    counts=pd.read_csv(counts_path,sep="\t",index_col=0)
    counts.index=counts.index.astype(str)
    pairs=pd.read_csv(pairs_path,sep="\t",dtype=str)
    samples=pd.read_csv(samples_path,sep="\t",dtype=str)
    for frame,columns,name in ((pairs,{"pair_id","gene1","gene2"},"pairs"),
                              (samples,{"sample_id","tissue","replicate"},"samples")):
        if not columns.issubset(frame): raise ValueError(f"Missing {name} columns: {columns-set(frame)}")
        if frame[list(columns)].isna().any().any(): raise ValueError(f"Missing {name} identifiers")
        if frame[list(columns)].apply(lambda x:x.str.strip().eq("")).any().any(): raise ValueError(f"Empty {name} identifiers")
    if pairs.pair_id.duplicated().any() or samples.sample_id.duplicated().any(): raise ValueError("Duplicate identifiers")
    if pairs.gene1.eq(pairs.gene2).any():raise ValueError("Self-pairs are invalid")
    if samples.duplicated(["tissue","replicate"]).any():raise ValueError("Replicate identifiers must be unique within tissue")
    if not samples.tissue.map(lambda s:bool(re.fullmatch(r"[A-Za-z0-9_-]+",s))).all():
        raise ValueError("Tissue identifiers must use letters, digits, underscore or hyphen")
    def check(frame,name):
        if frame.index.has_duplicates or frame.columns.has_duplicates:raise ValueError(f"Duplicate {name} rows or columns")
        if list(frame.columns)!=samples.sample_id.tolist():raise ValueError(f"{name} columns must match ordered sample metadata")
        values=frame.to_numpy(dtype=float)
        if not np.isfinite(values).all() or (values<0).any() or (values!=np.floor(values)).any():
            raise ValueError(f"{name} requires finite nonnegative integer counts, not TPM or fractional estimates")
        if (values>np.iinfo(np.int32).max).any():raise ValueError("Counts exceed DESeq2 integer capacity")
    check(counts,"whole-gene counts")
    missing=(set(pairs.gene1)|set(pairs.gene2))-set(counts.index)
    if missing:raise ValueError(f"Whole-gene counts lack {len(missing)} pair genes")
    if bool(copy1_path)!=bool(copy2_path):raise ValueError("Supply both copy count matrices or neither")
    if copy1_path:
        for path in (copy1_path,copy2_path):
            table=pd.read_csv(path,sep="\t",index_col=0);check(table,"pair-copy counts")
            if list(table.index.astype(str))!=pairs.pair_id.tolist():raise ValueError("Copy matrix pair rows must match pairs in order")
    if (samples.groupby("tissue").size()<2).any():raise ValueError("At least two biological libraries per tissue are required")
    return dict(genes=len(counts),pairs=len(pairs),libraries=len(samples),tissues=int(samples.tissue.nunique()))


def fit_copy_counts(counts_path,pairs_path,samples_path,output_dir,*,rscript="Rscript",test="wald",
                    copy1_path=None,copy2_path=None,observability_method="uncorrected",alpha=.05):
    if test not in {"wald","lrt"}:raise ValueError("Unsupported test")
    if observability_method not in {"reciprocal-common-cds","equal-effective-length","uncorrected"}:
        raise ValueError("Unknown observability method")
    if not 0<alpha<1:raise ValueError("Alpha must lie strictly between zero and one")
    if observability_method=="reciprocal-common-cds" and not (copy1_path and copy2_path):
        raise ValueError("Reciprocal mode requires the two pair-specific raw-count matrices")
    dimensions=validate_count_inputs(counts_path,pairs_path,samples_path,copy1_path,copy2_path)
    out=Path(output_dir).resolve()
    if out.exists():raise FileExistsError("Count fitting output must be new; previous results are never overwritten")
    out.mkdir(parents=True)
    files=dict(counts=counts_path,pairs=pairs_path,samples=samples_path,copy1=copy1_path,copy2=copy2_path)
    inputs={k:dict(path=str(Path(v).resolve()),sha256=hashlib.sha256(Path(v).read_bytes()).hexdigest()) for k,v in files.items() if v}
    request=dict(inputs=inputs,test=test,alpha=alpha,dimensions=dimensions,observability_method=observability_method,
        model="~biological_library+copy",normalization="Whole-gene poscounts size factors duplicated across copies",
        count_filter="sum across copies and libraries >=10",independent_filtering=False,
        multiple_testing="BH across tested pairs separately within each tissue; not final-class FDR",
        status="RUNNING",warning="Experimental count engine: high-variance few-replicate error control is not universally established; expression classes are not proven functional fates")
    request_path=out/"count_fit_metadata.json";request_path.write_text(json.dumps(request,indent=2))
    command=[rscript,str(Path(__file__).parent/"resources/fit_copy_counts.R"),
        str(Path(counts_path).resolve()),str(Path(pairs_path).resolve()),str(Path(samples_path).resolve()),str(out),test,
        str(Path(copy1_path).resolve()) if copy1_path else "NONE",str(Path(copy2_path).resolve()) if copy2_path else "NONE",str(alpha)]
    env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1",MKL_NUM_THREADS="1")
    try:
        with (out/"count_fit.log").open("w") as log:
            subprocess.run(command,stdout=log,stderr=log,check=True,env=env)
    except Exception as exc:
        request["status"]="FAILED";request["error"]=repr(exc)
        request_path.write_text(json.dumps(request,indent=2));raise
    request["status"]="COMPLETE";request_path.write_text(json.dumps(request,indent=2))
    return request_path
