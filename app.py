"""Streamlit interface for iDGRM's explicit two-stage workflow."""

from __future__ import annotations

import io
import tempfile
import zipfile
from pathlib import Path

import streamlit as st

from idgrm.config import IDGRMConfig
from idgrm.pipeline import run_refinement_analysis, run_science_analysis

st.set_page_config(page_title="iDGRM", page_icon="🧬", layout="wide")
st.title("🧬 iDGRM")
st.caption("先复现 Science 四分类，再对 AED 与 SUB_OR_NEO 做可追溯的扩展细分")
st.info(
    "第一阶段只输出 UNMAPPED、NO_DIFFERENCE、AED、SUB_OR_NEO。"
    "第二阶段不改变主分类，只为 AED 和 SUB_OR_NEO 增加二级标签。"
)


def _save(upload: object, folder: Path, fallback: str) -> Path | None:
    if upload is None:
        return None
    target = folder / fallback
    target.write_bytes(upload.getvalue())
    return target


def _archive(folder: Path) -> bytes:
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as archive:
        for path in folder.iterdir():
            if path.is_file():
                archive.write(path, arcname=path.name)
    return buffer.getvalue()


def _config(prefix: str) -> tuple[IDGRMConfig, str]:
    left, right = st.columns(2)
    with left:
        log2fc = st.number_input("|log2FC| 阈值", 0.01, value=1.0, step=0.1, key=f"{prefix}_lfc")
        min_expression = st.number_input(
            "表达 on/off 阈值", 0.000001, value=1.0, step=0.1, key=f"{prefix}_expr"
        )
        min_tissues = st.number_input("最少可评价组织数", 1, value=2, step=1, key=f"{prefix}_nt")
    with right:
        alpha = st.number_input(
            "P/q 阈值", 0.000001, 1.0, value=0.001, format="%.6f", key=f"{prefix}_alpha"
        )
        p_adjust = st.selectbox("多重检验", ["bh", "none"], key=f"{prefix}_padj")
        normalization = st.selectbox(
            "列归一化", ["none", "cpm", "median-ratio"], key=f"{prefix}_norm"
        )
    return IDGRMConfig(
        log2fc_threshold=float(log2fc), alpha=float(alpha), p_adjust=p_adjust,
        min_expression=float(min_expression), min_evaluable_tissues=int(min_tissues),
    ), normalization


science_tab, refine_tab = st.tabs(["第一阶段：Science 四分类", "第二阶段：iDGRM 扩展细分"])

with science_tab:
    st.markdown("输出严格限定为四类：**UNMAPPED、NO_DIFFERENCE、AED、SUB_OR_NEO**。")
    expression_file = st.file_uploader("表达矩阵", type=["csv", "tsv", "txt"], key="s_expr")
    pairs_file = st.file_uploader("重复基因对", type=["csv", "tsv", "txt"], key="s_pairs")
    samples_file = st.file_uploader("样本信息（可选）", type=["csv", "tsv", "txt"], key="s_samples")
    science_config, science_norm = _config("science")
    if st.button("运行 Science 四分类", type="primary", use_container_width=True):
        if expression_file is None or pairs_file is None:
            st.error("请上传表达矩阵和重复基因对。")
        else:
            try:
                with tempfile.TemporaryDirectory(prefix="idgrm-science-") as temp_name:
                    temp, output = Path(temp_name), Path(temp_name) / "science"
                    result = run_science_analysis(
                        expression_path=_save(expression_file, temp, "expression.tsv"),
                        pairs_path=_save(pairs_file, temp, "pairs.tsv"),
                        samples_path=_save(samples_file, temp, "samples.tsv"),
                        output_dir=output, config=science_config, normalization=science_norm,
                    )
                    st.success(f"完成：{len(result.classifications)} 对；主分类仅包含四类。")
                    st.dataframe(result.summary, use_container_width=True, hide_index=True)
                    st.dataframe(result.classifications[[
                        "pair_id", "gene1", "gene2", "science_class"
                    ]], use_container_width=True, hide_index=True)
                    st.download_button(
                        "下载第一阶段完整结果 ZIP", _archive(output),
                        "idgrm_science_four_class.zip", "application/zip", use_container_width=True,
                    )
            except Exception as exc:
                st.exception(exc)

with refine_tab:
    st.markdown(
        "上传第一阶段结果。AED细分为普通AED或表达丢失样AED；"
        "SUB_OR_NEO细分为SUB、NEO-like或歧义型。"
    )
    classification_file = st.file_uploader(
        "science_classifications.tsv", type=["csv", "tsv", "txt"], key="r_class"
    )
    evidence_file = st.file_uploader(
        "tissue_evidence.tsv", type=["csv", "tsv", "txt"], key="r_ev"
    )
    refine_pairs_file = st.file_uploader(
        "重复基因对（与第一阶段相同）", type=["csv", "tsv", "txt"], key="r_pairs"
    )
    ancestor_expression_file = st.file_uploader(
        "外群/祖先表达矩阵（推荐）", type=["csv", "tsv", "txt"], key="r_ancestor"
    )
    ancestor_samples_file = st.file_uploader(
        "外群样本信息（可选）", type=["csv", "tsv", "txt"], key="r_ancestor_samples"
    )
    refine_config, refine_norm = _config("refine")
    loss_fraction = st.number_input(
        "表达丢失样 AED 的组织比例阈值", 0.0, 1.0, value=0.8, step=0.05
    )
    refine_config = IDGRMConfig(**{
        **refine_config.to_dict(), "loss_fraction": float(loss_fraction)
    })
    if st.button("运行扩展细分", type="primary", use_container_width=True):
        if any(item is None for item in (classification_file, evidence_file, refine_pairs_file)):
            st.error("请上传第一阶段分类、组织证据和相同的重复基因对。")
        elif ancestor_samples_file is not None and ancestor_expression_file is None:
            st.error("上传外群样本信息时也必须上传外群表达矩阵。")
        else:
            try:
                with tempfile.TemporaryDirectory(prefix="idgrm-refine-") as temp_name:
                    temp, output = Path(temp_name), Path(temp_name) / "refined"
                    result = run_refinement_analysis(
                        science_classifications_path=_save(classification_file, temp, "science.tsv"),
                        evidence_path=_save(evidence_file, temp, "evidence.tsv"),
                        pairs_path=_save(refine_pairs_file, temp, "pairs.tsv"),
                        ancestor_expression_path=_save(ancestor_expression_file, temp, "ancestor.tsv"),
                        ancestor_samples_path=_save(ancestor_samples_file, temp, "ancestor_samples.tsv"),
                        output_dir=output, config=refine_config, normalization=refine_norm,
                    )
                    st.success(f"完成：细分 {len(result.classifications)} 个候选。")
                    st.dataframe(result.summary, use_container_width=True, hide_index=True)
                    st.dataframe(result.classifications[[
                        "pair_id", "parent_science_class", "extended_subtype",
                        "lost_copy", "innovating_copy", "refinement_reason",
                    ]], use_container_width=True, hide_index=True)
                    st.download_button(
                        "下载第二阶段完整结果 ZIP", _archive(output),
                        "idgrm_extended_refinement.zip", "application/zip", use_container_width=True,
                    )
                    if result.warnings:
                        st.warning("\n".join(result.warnings))
            except Exception as exc:
                st.exception(exc)
