"""Self-contained HTML reporting for iDGRM runs."""

from __future__ import annotations

import html
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd


COLORS = {
    "SUBFUNCTIONALIZATION": "#2a9d8f",
    "NEOFUNCTIONALIZATION": "#e76f51",
    "AMBIGUOUS_SUB_NEO": "#8d6cab",
    "AED": "#457b9d",
    "EXPRESSION_LOSS": "#6c757d",
    "NO_DIFFERENCE": "#b7b7a4",
    "INSUFFICIENT_DATA": "#d9d9d9",
}


def _value(value: object) -> str:
    if pd.isna(value):
        return ""
    if isinstance(value, float):
        return f"{value:.4g}"
    return str(value)


def write_html_report(
    path: str | Path,
    classifications: pd.DataFrame,
    summary: pd.DataFrame,
    evidence: pd.DataFrame,
    metadata: dict[str, Any],
) -> Path:
    """Write a portable report with no remote JavaScript or CSS dependencies."""

    path = Path(path)
    extended = summary[summary["scheme"].eq("extended_fate")].copy()
    total = max(1, int(len(classifications)))
    cards: list[str] = []
    bars: list[str] = []
    for row in extended.itertuples(index=False):
        color = COLORS.get(str(row.class_code), "#777")
        cards.append(
            "<div class='card'><span class='dot' style='background:%s'></span>"
            "<div><strong>%s</strong><small>%s</small></div><b>%s</b></div>"
            % (
                color,
                html.escape(str(row.label_cn)),
                html.escape(str(row.class_code)),
                html.escape(str(row.n_pairs)),
            )
        )
        width = max(0.5, float(row.proportion) * 100.0)
        bars.append(
            "<div class='bar-row'><span>%s</span><div class='track'>"
            "<div class='fill' style='width:%.3f%%;background:%s'></div></div>"
            "<b>%.1f%%</b></div>"
            % (html.escape(str(row.label_cn)), width, color, float(row.proportion) * 100.0)
        )

    display_columns = [
        "pair_id",
        "gene1",
        "gene2",
        "duplication_type",
        "science_class",
        "extended_fate",
        "evidence_basis",
        "confidence",
        "major_copy",
        "innovating_copy",
        "lost_copy",
        "n_gene1_high",
        "n_gene2_high",
        "classification_reason",
    ]
    display_columns = [column for column in display_columns if column in classifications.columns]
    table_rows: list[str] = []
    for row in classifications[display_columns].head(10_000).itertuples(index=False, name=None):
        values = dict(zip(display_columns, row))
        fate = str(values.get("extended_fate", ""))
        cells = "".join(
            f"<td>{html.escape(_value(values[column]))}</td>" for column in display_columns
        )
        table_rows.append(
            f"<tr data-fate='{html.escape(fate)}' style='--fate:{COLORS.get(fate, '#777')}'>{cells}</tr>"
        )
    headers = "".join(f"<th>{html.escape(column)}</th>" for column in display_columns)
    config_json = html.escape(
        json.dumps(metadata.get("config", {}), ensure_ascii=False, indent=2)
    )
    warning_items = "".join(
        f"<li>{html.escape(str(item))}</li>" for item in metadata.get("warnings", [])
    ) or "<li>None</li>"
    proxy_count = int(classifications["evidence_basis"].eq("expression_only_proxy").sum())
    generated = datetime.now(timezone.utc).isoformat(timespec="seconds")

    document = f"""<!doctype html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>iDGRM analysis report</title>
<style>
:root {{ --ink:#183153; --muted:#667085; --paper:#fff; --line:#e6eaf0; --bg:#f5f7fb; }}
* {{ box-sizing:border-box; }}
body {{ margin:0; background:var(--bg); color:var(--ink); font-family:Inter,Segoe UI,"Microsoft YaHei",sans-serif; }}
.wrap {{ width:min(1500px,96vw); margin:32px auto 64px; }}
.hero {{ background:linear-gradient(135deg,#173b57,#255f6f 58%,#2a9d8f); color:white; border-radius:20px; padding:30px 34px; box-shadow:0 16px 40px #173b5724; }}
.hero h1 {{ margin:0 0 8px; font-size:34px; letter-spacing:.02em; }}
.hero p {{ margin:5px 0; color:#e7f5f3; max-width:900px; }}
.hero .meta {{ font-size:13px; opacity:.86; margin-top:14px; }}
.section {{ background:var(--paper); border:1px solid var(--line); border-radius:18px; padding:24px; margin-top:20px; box-shadow:0 8px 26px #1831530a; }}
.section h2 {{ margin:0 0 18px; font-size:20px; }}
.cards {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(190px,1fr)); gap:12px; }}
.card {{ display:flex; gap:10px; align-items:center; border:1px solid var(--line); border-radius:13px; padding:13px; }}
.card .dot {{ width:12px; height:40px; border-radius:8px; flex:none; }}
.card div {{ min-width:0; flex:1; }} .card strong,.card small {{ display:block; }}
.card small {{ color:var(--muted); font-size:11px; margin-top:3px; }} .card b {{ font-size:22px; }}
.bars {{ margin-top:20px; display:grid; gap:10px; }}
.bar-row {{ display:grid; grid-template-columns:170px 1fr 60px; align-items:center; gap:10px; font-size:13px; }}
.track {{ height:12px; background:#eef1f5; border-radius:10px; overflow:hidden; }} .fill {{ height:100%; border-radius:10px; }}
.controls {{ display:flex; flex-wrap:wrap; gap:10px; margin-bottom:12px; }}
input,select {{ border:1px solid #cdd5df; border-radius:9px; padding:9px 11px; color:var(--ink); background:white; }}
input {{ min-width:320px; }}
.table-wrap {{ overflow:auto; max-height:68vh; border:1px solid var(--line); border-radius:12px; }}
table {{ border-collapse:separate; border-spacing:0; width:100%; font-size:12px; }}
th {{ position:sticky; top:0; z-index:1; background:#f0f4f8; text-align:left; padding:10px; border-bottom:1px solid var(--line); white-space:nowrap; }}
td {{ padding:9px 10px; border-bottom:1px solid #edf0f4; vertical-align:top; max-width:430px; }}
tbody tr {{ border-left:4px solid var(--fate); }} tbody tr:hover {{ background:#f8fbfd; }}
pre {{ background:#101828; color:#d8e4ed; border-radius:12px; padding:16px; overflow:auto; }}
.note {{ color:var(--muted); line-height:1.65; }}
.pill {{ display:inline-block; border-radius:999px; background:#e9f5f3; color:#176b61; padding:5px 10px; font-size:12px; margin-right:6px; }}
ul {{ line-height:1.65; }}
@media (max-width:700px) {{ .bar-row {{ grid-template-columns:110px 1fr 48px; }} input {{ min-width:100%; }} }}
</style>
</head>
<body><main class="wrap">
<section class="hero">
  <h1>iDGRM</h1>
  <p>Inference of Duplicated-Gene Retention Mechanisms · 重复基因表达命运分类报告</p>
  <p class="meta">Generated {html.escape(generated)} · {total} gene pairs · {len(evidence)} tissue-evidence rows</p>
</section>
<section class="section">
  <h2>分类概览</h2>
  <div class="cards">{''.join(cards)}</div>
  <div class="bars">{''.join(bars)}</div>
</section>
<section class="section">
  <h2>证据解释</h2>
  <p class="note"><span class="pill">Science-compatible</span>science_class 保留论文的 reciprocal sub/neo、AED 和 no-difference 逻辑。</p>
  <p class="note"><span class="pill">Extended fate</span>extended_fate 将 reciprocal 候选拆分为 sub、neo 或 ambiguous，并增加表达丢失与数据不足类别。</p>
  <p class="note">本次有 <strong>{proxy_count}</strong> 对使用 expression-only proxy。没有外群单拷贝正交基因时，“neo”表示受限的新表达域代理，不能单凭表达证明蛋白质层面的新功能。</p>
</section>
<section class="section">
  <h2>基因对结果</h2>
  <div class="controls">
    <input id="search" type="search" placeholder="搜索 pair / gene / reason..." oninput="filterRows()">
    <select id="fate" onchange="filterRows()"><option value="">全部类别</option>{''.join(f'<option value="{html.escape(code)}">{html.escape(code)}</option>' for code in COLORS)}</select>
  </div>
  <div class="table-wrap"><table><thead><tr>{headers}</tr></thead><tbody id="rows">{''.join(table_rows)}</tbody></table></div>
  <p class="note">HTML 最多展示 10,000 对；完整结果见 classifications.tsv。</p>
</section>
<section class="section">
  <h2>运行参数</h2><pre>{config_json}</pre>
  <h2>警告</h2><ul>{warning_items}</ul>
</section>
</main>
<script>
function filterRows() {{
  const query = document.getElementById('search').value.toLowerCase();
  const fate = document.getElementById('fate').value;
  document.querySelectorAll('#rows tr').forEach(row => {{
    const matchText = !query || row.innerText.toLowerCase().includes(query);
    const matchFate = !fate || row.dataset.fate === fate;
    row.style.display = matchText && matchFate ? '' : 'none';
  }});
}}
</script></body></html>"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(document, encoding="utf-8")
    return path
