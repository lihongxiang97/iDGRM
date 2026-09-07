"""Streamlit browser interface for iDGRM."""

from __future__ import annotations

import io
import tempfile
import zipfile
from pathlib import Path

import streamlit as st

from idgrm.config import IDGRMConfig
from idgrm.pipeline import run_analysis


st.set_page_config(page_title="iDGRM", page_icon="🧬", layout="wide")
st.title("🧬 iDGRM")
st.caption("Inference of Duplicated-Gene Retention Mechanisms · 基于多组织表达的重复基因进化命运分类")

with st.expander("输入格式与判定边界", expanded=False):
    st.markdown(
        """
- 表达矩阵：第一列为 `gene_id`，其余列为样本；若不上传样本信息，则其余列按组织均值处理。
- 样本信息：至少包含 `sample_id` 与 `tissue`。
- 基因对：支持 `gene1/gene2`、`Duplicate 1/2` 或 `Parent copy/Daughter copy`。
- 若提供外群单拷贝正交基因表达及 `ancestor_id`，sub/neo 使用祖先参照；否则结果明确标记为 expression-only proxy。
        """
    )

left, right = st.columns([1.15, 0.85])
with left:
    expression_file = st.file_uploader("表达矩阵（CSV/TSV）", type=["csv", "tsv", "txt"])
    pairs_file = st.file_uploader("重复基因对（CSV/TSV）", type=["csv", "tsv", "txt"])
    samples_file = st.file_uploader("样本信息（可选）", type=["csv", "tsv", "txt"])
    ancestor_expression_file = st.file_uploader("外群/祖先表达矩阵（可选）", type=["csv", "tsv", "txt"])
    ancestor_samples_file = st.file_uploader("外群样本信息（可选）", type=["csv", "tsv", "txt"])

with right:
    st.subheader("参数")
    log2fc = st.number_input("|log2FC| 阈值", min_value=0.01, value=1.0, step=0.1)
    alpha = st.number_input("P/q 阈值", min_value=0.000001, max_value=1.0, value=0.001, format="%.6f")
    min_expression = st.number_input(
        "表达 on/off 阈值", min_value=0.000001, value=1.0, step=0.1, format="%.6f"
    )
    min_tissues = st.number_input("最少可评价组织数", min_value=1, value=2, step=1)
    p_adjust = st.selectbox("多重检验", ["bh", "none"], index=0)
    normalization = st.selectbox("列归一化", ["none", "cpm", "median-ratio"], index=0)
    detect_loss = st.checkbox("识别表达丢失", value=True)


def _save(upload: object, folder: Path, fallback: str) -> Path | None:
    if upload is None:
        return None
    # Fixed role-specific names prevent two uploads with the same local filename
    # from silently overwriting each other in the temporary run directory.
    target = folder / fallback
    target.write_bytes(upload.getvalue())
    return target


if st.button("运行 iDGRM", type="primary", use_container_width=True):
    if expression_file is None or pairs_file is None:
        st.error("请至少上传表达矩阵和重复基因对。")
    elif ancestor_samples_file is not None and ancestor_expression_file is None:
        st.error("上传外群样本信息时也必须上传外群表达矩阵。")
    else:
        try:
            with st.spinner("正在构建组织证据并分类..."):
                with tempfile.TemporaryDirectory(prefix="idgrm-") as temp_name:
                    temp = Path(temp_name)
                    expression_path = _save(expression_file, temp, "expression.tsv")
                    pairs_path = _save(pairs_file, temp, "pairs.tsv")
                    samples_path = _save(samples_file, temp, "samples.tsv")
                    ancestor_expression_path = _save(
                        ancestor_expression_file, temp, "ancestor_expression.tsv"
                    )
                    ancestor_samples_path = _save(
                        ancestor_samples_file, temp, "ancestor_samples.tsv"
                    )
                    output = temp / "results"
                    config = IDGRMConfig(
                        log2fc_threshold=float(log2fc),
                        alpha=float(alpha),
                        p_adjust=p_adjust,
                        min_expression=float(min_expression),
                        min_evaluable_tissues=int(min_tissues),
                        detect_expression_loss=detect_loss,
                    )
                    result = run_analysis(
                        expression_path=expression_path,
                        pairs_path=pairs_path,
                        samples_path=samples_path,
                        ancestor_expression_path=ancestor_expression_path,
                        ancestor_samples_path=ancestor_samples_path,
                        output_dir=output,
                        config=config,
                        normalization=normalization,
                    )
                    st.success(f"完成：{len(result.classifications)} 对重复基因")
                    st.subheader("分类概览")
                    st.dataframe(
                        result.summary[result.summary["scheme"].eq("extended_fate")],
                        use_container_width=True,
                        hide_index=True,
                    )
                    st.subheader("结果预览")
                    preview_columns = [
                        "pair_id", "gene1", "gene2", "science_class", "extended_fate",
                        "evidence_basis", "confidence", "classification_reason",
                    ]
                    st.dataframe(
                        result.classifications[preview_columns],
                        use_container_width=True,
                        hide_index=True,
                    )
                    archive_buffer = io.BytesIO()
                    with zipfile.ZipFile(archive_buffer, "w", zipfile.ZIP_DEFLATED) as archive:
                        for file_path in output.iterdir():
                            archive.write(file_path, arcname=file_path.name)
                    st.download_button(
                        "下载完整结果 ZIP",
                        data=archive_buffer.getvalue(),
                        file_name="idgrm_results.zip",
                        mime="application/zip",
                        use_container_width=True,
                    )
                    if result.warnings:
                        st.warning("\n".join(result.warnings))
        except Exception as exc:
            st.exception(exc)
