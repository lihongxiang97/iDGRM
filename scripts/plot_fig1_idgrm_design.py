"""Draw the iDGRM design schematic using editable Matplotlib vectors.

All profiles are illustrative, not experimental observations.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["font.size"] = 7

import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle, Circle, PathPatch
from matplotlib.path import Path as MPath

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from idgrm.classifier import classify_primary_pairs
from idgrm.config import IDGRMConfig
from idgrm import __version__

W, H = 183., 180.
C = dict(ink="#243444", muted="#647483", line="#D7DFE5", pale="#F5F8FA",
         blue="#3E7396", rose="#B27687", violet="#78658F", teal="#4F8C8A",
         gray="#939EA7", lightblue="#EDF3F7", lightviolet="#F3F0F7")
TEXT = []


def text(ax, x, y, label, size=7, color=None, weight="normal", ha="left", va="top", **kwargs):
    artist = ax.text(x, y, label, fontsize=size, color=color or C["ink"],
        fontweight=weight, ha=ha, va=va, linespacing=1.35, **kwargs)
    TEXT.append(artist)
    return artist


def line(ax, xs, ys, color=None, lw=.7, **kwargs):
    ax.plot(xs, ys, color=color or C["line"], lw=lw, solid_capstyle="round", **kwargs)


def box(ax, x, y, w, h, face="white", edge=None, radius=1.5, lw=.7, **kwargs):
    patch = FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={radius}",
        facecolor=face, edgecolor=edge or C["line"], linewidth=lw, **kwargs)
    ax.add_patch(patch)
    return patch


def arrow(ax, start, end, color=None, dashed=False, lw=.9, scale=8):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=scale,
        linewidth=lw, color=color or C["muted"], shrinkA=0, shrinkB=0,
        linestyle=(0, (3, 2)) if dashed else "-"))


def panel(ax, letter, x, y, title, tag=None):
    text(ax, x, y-.2, letter, 10, weight="bold")
    text(ax, x+6, y, title, 8.5, weight="bold")
    if tag:
        text(ax, 175, y+.6, tag, 6.2, color=C["muted"], ha="right")


def gene(ax, x, y, color, scale=1):
    line(ax, [x, x+12*scale], [y, y], color, .9)
    for start, width in [(1, 2), (5, 2.2), (9, 2)]:
        ax.add_patch(Rectangle((x+start*scale, y-1.05*scale), width*scale, 2.1*scale,
            facecolor=color, edgecolor="none"))
    arrow(ax, (x+11*scale, y), (x+13*scale, y), color, scale=5)


def matrix_icon(ax, x, y):
    colors = ["#C8D8E4", "#91B0C8", "#DCE7EE", "#5C87A6"]
    for i in range(3):
        for j in range(5):
            ax.add_patch(Rectangle((x+j*1.9, y+i*1.9), 1.45, 1.45,
                facecolor=colors[(i*3+j) % 4], edgecolor="none"))


def profiles():
    return [
        dict(key="unmappable", title="Unmappable", color=C["gray"],
            a=[np.nan]*4, b=[np.nan]*4, direction=["no_difference"]*4,
            rule="Insufficient callable evidence\n(default: fewer than 2 tissues)", expected="UNMAPPED"),
        dict(key="no_difference", title="No difference", color=C["teal"],
            a=[4.2, 5.2, 3.6, 5.5], b=[4.5, 5.0, 3.8, 5.2], direction=["no_difference"]*4,
            rule="Neither directional\nclassification rule is met", expected="NO_DIFFERENCE"),
        dict(key="asymmetric", title="Asymmetrically\nexpressed", color=C["blue"],
            a=[8.0, 6.8, 8.5, 6.0], b=[1.8, 1.7, 2.0, 1.5], direction=["gene1_high"]*4,
            rule="Same-copy dominance in ≥⅓ of N;\nno significant reversal", expected="ASYMMETRICALLY_EXPRESSED"),
        dict(key="subneo", title="Sub-/neofunctionalized", color=C["violet"],
            a=[8.2, 7.6, 1.6, 1.7], b=[1.6, 1.8, 7.8, 8.1],
            direction=["gene1_high", "gene1_high", "gene2_high", "gene2_high"],
            rule="Each copy dominates\nat least one tissue", expected="SUB_OR_NEOFUNCTIONALIZED"),
    ]


def check_rules(profile_specs):
    pairs, evidence = [], []
    for item in profile_specs:
        pairs.append(dict(pair_id=item["key"], gene1="a_"+item["key"], gene2="b_"+item["key"],
            duplication_type="illustration", role1="", role2="", ancestor_id=""))
        for i, (a, b, direction) in enumerate(zip(item["a"], item["b"], item["direction"])):
            evidence.append(dict(pair_id=item["key"], tissue=f"T{i+1}", mean1=a, mean2=b,
                active1=bool(np.isfinite(a)), active2=bool(np.isfinite(b)), direction=direction,
                log2fc_gene2_over_gene1=np.log2(b/a) if np.isfinite(a) else np.nan,
                pvalue=.8 if direction=="no_difference" else .001,
                qvalue=.8 if direction=="no_difference" else .001,
                method="ILLUSTRATIVE_RULE_INPUT_NOT_MEASURED"))
    classified = classify_primary_pairs(pd.DataFrame(evidence), pd.DataFrame(pairs), IDGRMConfig())
    labels = classified.set_index("pair_id").primary_class.to_dict()
    assert all(labels[p["key"]] == p["expected"] for p in profile_specs)
    return labels


def architecture(ax):
    panel(ax, "a", 8, 29, "From expression evidence to reproducible classifications", "CORE WORKFLOW")
    # Inputs are deliberately separated from the software engine.
    box(ax, 8, 40, 43, 43)
    text(ax, 12, 43.4, "USER-PREPARED INPUTS", 6.1, color=C["muted"], weight="bold")
    gene(ax, 12, 53, C["blue"], .72)
    gene(ax, 12, 56.6, C["rose"], .72)
    text(ax, 25, 52, "Duplicate pairs", 7, weight="bold")
    matrix_icon(ax, 12, 62)
    text(ax, 25, 62.2, "Integer copy counts", 7, weight="bold")
    text(ax, 25, 66, "+ whole-gene counts", 6.2, color=C["muted"])
    for i, height in enumerate([2, 4, 3, 5]):
        ax.add_patch(Rectangle((12+i*2.4, 78-height), 1.5, height, facecolor=C["teal"], edgecolor="none"))
    text(ax, 25, 73.6, "TPM + tissue metadata", 6.8, weight="bold")
    text(ax, 25, 77.4, "Biological replicates", 6.1, color=C["muted"])
    arrow(ax, (52.5, 61.5), (61, 61.5), C["blue"], scale=9)

    box(ax, 63, 40, 63, 43, C["lightblue"], edge="#CDDBE5")
    text(ax, 67, 43.4, "iDGRM STATISTICAL ENGINE", 6.1, color=C["blue"], weight="bold")
    text(ax, 67, 50.1, "Paired negative-binomial model", 7.3, weight="bold")
    text(ax, 67, 55.4, "DESeq2 in R  ·  ~ library + copy", 7, color=C["muted"])
    line(ax, [67, 122], [61, 61], "#CEDDE7", .6)
    text(ax, 67, 64.1, "BH-adjusted P ≤ 0.05", 7.3, weight="bold")
    text(ax, 67, 69.1, "|log2FC| ≥ 1", 7.3, weight="bold")
    text(ax, 67, 76.3, "Activity + valid-test filtering", 6.8, color=C["muted"])
    arrow(ax, (127.5, 61.5), (134, 61.5), C["blue"], scale=9)

    output_specs = [
        ("Unmappable", C["gray"]), ("No difference", C["teal"]),
        ("Asymmetrically expressed", C["blue"]), ("Sub-/neofunctionalized", C["violet"])]
    for i, (label, color) in enumerate(output_specs):
        y = 41.3+i*10.6
        box(ax, 136, y, 39, 7.7, "white", edge="#D9E1E6", radius=1.1, lw=.65)
        ax.add_patch(Rectangle((136, y+1.5), .8, 4.7, facecolor=color, edgecolor="none"))
        text(ax, 139, y+2.25, label, 6.7, color=color, weight="bold")

    text(ax, 8, 87, "Input prerequisite", 6.5, weight="bold", color=C["muted"])
    text(ax, 34, 87, "Copy-comparable quantification; strict, uniquely assignable reads", 6.5, color=C["muted"])
    line(ax, [8, 175], [94, 94], C["line"], .65)
    text(ax, 8, 96.5, "Flexible entry", 6.3, weight="bold", color=C["blue"])
    text(ax, 31.5, 96.5, "Import completed DE tests, or use the separate log-expression paired t-test route", 6.3, color=C["muted"])


def primary_logic(ax, specs):
    panel(ax, "b", 8, 108, "Cross-tissue decision rules")
    text(ax, 5.1, 141.5, "Expression (a.u.)", 5.9, color=C["muted"], rotation=90, ha="center", va="center")
    # One shared copy identity key; line style and marker shape retain grayscale meaning.
    line(ax, [122, 128], [110, 110], C["blue"], 1.2)
    ax.plot([125], [110], "o", ms=2.7, color=C["blue"])
    text(ax, 130, 108.8, "Copy 1", 6.5)
    line(ax, [148, 154], [110, 110], C["rose"], 1.2, linestyle=(0, (2, 1)))
    ax.plot([151], [110], "s", ms=2.5, color=C["rose"])
    text(ax, 156, 108.8, "Copy 2", 6.5)
    for i, item in enumerate(specs):
        x = 8+i*43
        line(ax, [x, x+38], [118.7, 118.7], item["color"], 1.5)
        text(ax, x, 121.2, item["title"], 7.0, color=item["color"], weight="bold")
        px, py, pw, ph = x+2, 133, 32, 17
        line(ax, [px, px, px+pw], [py, py+ph, py+ph], C["line"], .65)
        xx = np.linspace(px+2, px+pw-1, 4)
        if i == 0:
            for j in range(4):
                box(ax, xx[j]-1.6, py+5.7, 3.2, 6, C["pale"], edge=C["line"], radius=.4, hatch="///", lw=.4)
            text(ax, px+pw/2, py+1, "Not evaluable", 6.5, color=C["muted"], ha="center")
        else:
            for values, color, marker, style in [(item["a"], C["blue"], "o", "-"),
                    (item["b"], C["rose"], "s", (0, (2, 1)))]:
                yy = py+ph-np.array(values)/10*ph
                line(ax, xx, yy, color, 1.2, linestyle=style)
                ax.plot(xx, yy, linestyle="none", marker=marker, ms=2.6, color=color)
        for j, xj in enumerate(xx):
            text(ax, xj, 151.6, f"T{j+1}", 5.9, color=C["muted"], ha="center")
        text(ax, x, 157.2, item["rule"], 6.4, color=C["ink"])
    text(ax, 8, 168.8, "Illustrative profiles, not measurements. N = evaluable active tissues. No difference does not imply functional equivalence.",
        6.15, color=C["muted"])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="figures/Fig1_iDGRM_design")
    args = parser.parse_args()
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    specs = profiles()
    checked = check_rules(specs)
    fig, ax = plt.subplots(figsize=(W/25.4, H/25.4), dpi=150)
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    ax.set(xlim=(0, W), ylim=(H, 0))
    ax.set_axis_off()
    text(ax, 8, 7, "iDGRM", 22, weight="bold", color=C["ink"])
    text(ax, 175, 11, "DESIGN  /  PRINCIPLES  /  WORKFLOW", 6.5, color=C["muted"], ha="right")
    text(ax, 8, 19, "Four-class classification of duplicate-gene expression patterns across tissues", 7.6, color=C["muted"])
    architecture(ax)
    primary_logic(ax, specs)
    text(ax, 8, 175, f"iDGRM {__version__} · Expression-pattern classes; SUB and NEO remain combined.", 6.2, color=C["muted"])
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    bounds = fig.bbox
    outside = [t.get_text() for t in TEXT if not bounds.contains(*t.get_window_extent(renderer).get_points()[0])
        or not bounds.contains(*t.get_window_extent(renderer).get_points()[1])]
    collisions = []
    for i, first in enumerate(TEXT):
        a = first.get_window_extent(renderer)
        for second in TEXT[i+1:]:
            b = second.get_window_extent(renderer)
            overlap_w = min(a.x1, b.x1)-max(a.x0, b.x0)
            overlap_h = min(a.y1, b.y1)-max(a.y0, b.y0)
            if overlap_w > 1 and overlap_h > 1:
                collisions.append([first.get_text(), second.get_text()])
    base = out / "Fig1_iDGRM_design"
    fig.savefig(base.with_suffix(".svg"), facecolor="white")
    svg_path = base.with_suffix(".svg")
    svg_path.write_text("\n".join(row.rstrip() for row in svg_path.read_text(encoding="utf-8").splitlines()) + "\n", encoding="utf-8")
    fig.savefig(base.with_suffix(".pdf"), facecolor="white")
    fig.savefig(base.with_suffix(".png"), dpi=300, facecolor="white")
    fig.savefig(out / "Fig1_iDGRM_design_preview.png", dpi=180, facecolor="white")
    fig.savefig(base.with_suffix(".tiff"), dpi=600, facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    svg_root = ET.parse(base.with_suffix(".svg")).getroot()
    n_svg_text = len(list(svg_root.iter("{http://www.w3.org/2000/svg}text")))
    with (out / "illustrative_profiles.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["illustration", "tissue", "copy1_arbitrary_units", "copy2_arbitrary_units", "illustrative_direction", "data_type"])
        for item in specs:
            for j, (a, b, d) in enumerate(zip(item["a"], item["b"], item["direction"])):
                writer.writerow([item["key"], f"T{j+1}", a, b, d, "ILLUSTRATIVE_NOT_MEASURED"])
    qa = dict(software_version=__version__, dimensions_mm=[W, H], backend="Python/Matplotlib",
        svg_editable_text_count=n_svg_text, out_of_canvas_text=outside, text_collisions=collisions,
        checked_primary_labels=checked, raster_dpi={"png": 300, "tiff": 600},
        integrity_note="All profiles are explicitly illustrative, not experimental data or accuracy evidence")
    qa["font"] = "Arial with sans-serif fallbacks; editable SVG text and embedded TrueType PDF"
    qa["qa_checks"] = {"in_canvas_text": not outside, "nonoverlapping_text": not collisions,
        "editable_svg_text": n_svg_text > 0, "code_checked_schematic_patterns": True}
    qa["sha256"] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir() if p.suffix in {".svg", ".pdf", ".png", ".tiff"}}
    (out / "figure_qa.json").write_text(json.dumps(qa, indent=2), encoding="utf-8")
    assert all(qa["qa_checks"].values()), "Figure QA failed; inspect figure_qa.json"
    print(json.dumps(qa, indent=2))


if __name__ == "__main__":
    main()
