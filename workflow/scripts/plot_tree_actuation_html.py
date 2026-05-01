#!/usr/bin/env python3
"""Interactive phylogenetic tree + FIRE actuation + MSA visualization via Plotly.

Layout (left to right):
  tree | FIRE % actuation | MSA (toggleable) with sequence logo above

Always writes HTML. Also writes PDF via kaleido (already in pixi deps).
The MSA panel can be shown/hidden with a button in the HTML toolbar.

Usage:
    python plot_tree_actuation_html.py \\
        --tree      trees/chr8_8453862_8454062.treefile \\
        --actuation fastas/chr8_8453862_8454062.tsv \\
        --alignment aligned/chr8_8453862_8454062.fa \\
        --html      plots/chr8_8453862_8454062.html \\
        --pdf       plots/chr8_8453862_8454062.pdf
"""

import argparse
import csv
import sys
from pathlib import Path

from Bio import Phylo
import plotly.graph_objects as go

# ── Colour palettes ───────────────────────────────────────────────────────────

# ColorBrewer Reds 10-bin: white (0%) → dark red (100%)
ACT_COLORS = [
    "#fff5f0", "#fee0d2", "#fcbba1", "#fc9272", "#fb6a4a",
    "#ef3b2c", "#cb181d", "#a50f15", "#67000d", "#3d0005",
]
ACT_LABELS = [f"{i*10}–{i*10+10}%" for i in range(10)]
ACT_MISSING = "#d0d0d0"   # FIRE actuation missing (light grey)

# DNA nucleotide colours (ClustalX-inspired) — only used for variable positions.
DNA_COLORS = {
    "A": "#009600",   # green
    "T": "#FF4444",   # red
    "G": "#FFB300",   # amber
    "C": "#4444FF",   # blue
}
MSA_INVARIANT = "#DEDEDE"   # light grey — base matches column consensus
MSA_MISSING   = "#888888"   # dark grey  — gap / N / unmapped

# MSA z-encoding:  0=var-A  1=var-T  2=var-G  3=var-C  4=invariant  5=missing
# Colorscale: one entry per integer z value, normalized to [0,1] via zmax=5.
# Must start at exactly 0.0 and end at exactly 1.0 — Plotly silently rejects
# colorscales that don't span [0,1] and falls back to "Reds".
BASE_TO_CODE = {"A": 0, "T": 1, "G": 2, "C": 3}
DNA_COLORSCALE = [
    [0.0, "#009600"],   # z=0: variant A — green
    [0.2, "#FF4444"],   # z=1: variant T — red
    [0.4, "#FFB300"],   # z=2: variant G — amber
    [0.6, "#4444FF"],   # z=3: variant C — blue
    [0.8, MSA_INVARIANT],   # z=4: matches consensus — light grey
    [1.0, MSA_MISSING],     # z=5: gap / N — dark grey
]


# ── Data loaders ──────────────────────────────────────────────────────────────

def load_actuation(path: Path) -> tuple[dict[str, float | None], dict[str, int | None]]:
    """Return (actuation_dict, coverage_dict) keyed by sample name."""
    actuation: dict[str, float | None] = {}
    coverage: dict[str, int | None] = {}
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            name = row["sample"]
            val = row.get("actuation", "NA")
            actuation[name] = None if val in ("NA", "") else float(val)
            cov = row.get("coverage", "NA")
            coverage[name] = None if cov in ("NA", "") else int(cov)
    return actuation, coverage


def load_alignment(path: Path) -> dict[str, str]:
    """Read a FASTA into {name: sequence}.

    Strips the leading ``_R_`` that MAFFT prepends when it reverse-complements
    a sequence (``--adjustdirection``), so names match the actuation TSV.
    """
    seqs: dict[str, str] = {}
    name = None
    parts: list[str] = []
    with path.open() as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(parts)
                raw = line[1:].split()[0]
                name = raw[3:] if raw.startswith("_R_") else raw
                parts = []
            else:
                parts.append(line)
    if name is not None:
        seqs[name] = "".join(parts)
    return seqs


# ── Tree position computation ─────────────────────────────────────────────────

def compute_positions(root) -> dict[int, dict]:
    """Assign x (cumulative branch length) and y (leaf index, 0-based) to every node."""
    pos: dict[int, dict] = {}

    def assign_y(clade, counter=[0]) -> float:
        if clade.is_terminal():
            y = float(counter[0])
            pos[id(clade)] = {"y": y, "name": clade.name, "terminal": True}
            counter[0] += 1
            return y
        child_ys = [assign_y(c) for c in clade.clades]
        y = (child_ys[0] + child_ys[-1]) / 2.0
        pos[id(clade)] = {"y": y, "terminal": False}
        return y

    def assign_x(clade, x: float = 0.0) -> None:
        pos[id(clade)]["x"] = x
        for child in clade.clades:
            assign_x(child, x + (child.branch_length or 0.0))

    assign_y(root)
    assign_x(root)
    return pos


def branch_xy(root, pos) -> tuple[list, list]:
    """Return x/y arrays (None-separated segments) for all tree branches."""
    xs: list = []
    ys: list = []
    for clade in root.find_clades(order="level"):
        if not clade.clades:
            continue
        px, py = pos[id(clade)]["x"], pos[id(clade)]["y"]
        child_ys = [pos[id(c)]["y"] for c in clade.clades]
        xs += [px, px, None]
        ys += [min(child_ys), max(child_ys), None]
        for child in clade.clades:
            xs += [px, pos[id(child)]["x"], None]
            ys += [pos[id(child)]["y"], pos[id(child)]["y"], None]
    return xs, ys


# ── Sequence logo helpers ─────────────────────────────────────────────────────

def logo_freqs(raw_z: list[list[int]], aln_len: int) -> dict[str, list[float]]:
    """Per-position frequency of A/T/G/C (gaps/N excluded from denominator).

    raw_z uses BASE_TO_CODE encoding (A=0,T=1,G=2,C=3; anything ≥4 = gap/missing).
    """
    freqs: dict[str, list[float]] = {nuc: [] for nuc in "ATGC"}
    for col in range(aln_len):
        non_gap = sum(1 for row in raw_z if row[col] < 4)
        for nuc, idx in BASE_TO_CODE.items():
            count = sum(1 for row in raw_z if row[col] == idx)
            freqs[nuc].append(count / non_gap if non_gap else 0.0)
    return freqs


# ── Main figure builder ───────────────────────────────────────────────────────

def make_figure(
    tree,
    actuation: dict[str, float | None],
    coverage: dict[str, int | None],
    alignment: dict[str, str],
    peak_id: str,
) -> go.Figure:
    tree.root_at_midpoint()
    tree.ladderize()

    # MAFFT --adjustdirection prefixes reversed-complemented leaf names with "_R_".
    # IQ-TREE inherits those names; strip the prefix so they match the actuation TSV.
    for _leaf in tree.get_terminals():
        if _leaf.name and _leaf.name.startswith("_R_"):
            _leaf.name = _leaf.name[3:]

    pos    = compute_positions(tree.root)
    leaves = list(tree.get_terminals())
    n      = len(leaves)
    names  = [lf.name for lf in leaves]
    ly     = [pos[id(lf)]["y"] for lf in leaves]   # 0 … n-1
    lx     = [pos[id(lf)]["x"] for lf in leaves]

    # Actuation
    acts   = [actuation.get(nm) for nm in names]
    covs   = [coverage.get(nm) for nm in names]
    colors = [
        ACT_COLORS[min(int(v * 10), 9)] if v is not None else ACT_MISSING
        for v in acts
    ]
    cov_labels = [f"#{c}" if c is not None else "" for c in covs]
    act_hover = [
        f"<b>{nm}</b><br>Actuation: {f'{v:.1%}' if v is not None else 'missing'}"
        f"<br>Coverage: {c if c is not None else 'missing'}"
        for nm, v, c in zip(names, acts, covs)
    ]

    # MSA — three-way encoding: variable ATGC | invariant | missing
    aln_len = max((len(s) for s in alignment.values()), default=0)
    raw_seqs = [alignment.get(nm, "-" * aln_len) for nm in names]

    # raw base codes: A/T/G/C → 0-3; everything else (gap/N) → 5
    raw_z = [
        [BASE_TO_CODE.get(c.upper(), 5) for c in seq]
        for seq in raw_seqs
    ]

    # Per-column consensus: majority non-gap base (code 0-3), or -1 if all gaps
    _code_to_char = ["A", "T", "G", "C"]
    consensus: list[int] = []
    for col in range(aln_len):
        counts = [0, 0, 0, 0]
        for row in raw_z:
            if row[col] < 4:
                counts[row[col]] += 1
        if any(c > 0 for c in counts):
            consensus.append(max(range(4), key=lambda i: counts[i]))  # noqa: B023
        else:
            consensus.append(-1)

    # Per-cell encoding: matches consensus → 4 (grey), differs → 0-3 (ATGC), gap → 5
    msa_z: list[list[int]] = []
    msa_cell: list[list[str]] = []   # single letter at variant positions only
    msa_hover: list[list[str]] = []  # rich hover tooltip
    for nm, seq, row_raw in zip(names, raw_seqs, raw_z):
        row_z, row_cell, row_hover = [], [], []
        for col_i, (c, code) in enumerate(zip(seq, row_raw)):
            cu = c.upper()
            cons = consensus[col_i]
            if code == 5:                       # gap / N / unmapped
                row_z.append(5)
                row_cell.append("")
                row_hover.append(f"<b>{nm}</b><br>pos {col_i+1}: {cu} (missing)")
            elif code == cons:                  # matches consensus — grey
                row_z.append(4)
                row_cell.append("")
                row_hover.append(f"<b>{nm}</b><br>pos {col_i+1}: {cu} (consensus)")
            else:                               # differs from consensus — ATGC colour
                row_z.append(code)
                row_cell.append(cu)
                cons_char = _code_to_char[cons] if cons >= 0 else "?"
                row_hover.append(
                    f"<b>{nm}</b><br>pos {col_i+1}: {cu} (variant, cons={cons_char})"
                )
        msa_z.append(row_z)
        msa_cell.append(row_cell)
        msa_hover.append(row_hover)

    freqs = logo_freqs(raw_z, aln_len)
    aln_positions = list(range(aln_len))

    # ── Axis domains ──────────────────────────────────────────────────────────
    # Horizontal: tree | labels | actuation | MSA
    # Vertical:   logo (top 8%) | main (bottom 89%)
    max_depth = max((pos[id(lf)]["x"] for lf in leaves), default=1.0)
    # Domains when MSA is shown (default)
    TREE_X  = [0.00, 0.25]
    LABEL_X = [0.26, 0.42]
    ACT_X   = [0.43, 0.50]
    MSA_X   = [0.52, 1.00]
    MAIN_Y  = [0.00, 0.89]
    LOGO_Y  = [0.92, 1.00]
    # Domains when MSA is hidden (tree+labels+actuation expand to fill)
    TREE_X_H  = [0.00, 0.38]
    LABEL_X_H = [0.39, 0.72]
    ACT_X_H   = [0.73, 0.88]

    # Shared y-range — explicit to eliminate auto-padding that breaks alignment
    Y_RANGE = [n - 0.5, -0.5]   # reversed: leaf 0 at top

    fig = go.Figure()

    # ── 1. Tree branches ──────────────────────────────────────────────────────
    bx, by = branch_xy(tree.root, pos)
    fig.add_trace(go.Scatter(
        x=bx, y=by,
        mode="lines",
        line=dict(color="black", width=0.9),
        hoverinfo="skip",
        showlegend=False,
        xaxis="x1", yaxis="y1",
    ))

    # ── 1b. Bootstrap support labels on internal nodes (≥70 only) ───────────
    boot_x, boot_y, boot_text = [], [], []
    for clade in tree.root.find_clades(order="level"):
        if not clade.is_terminal() and clade.confidence is not None:
            val = int(round(clade.confidence))
            if val >= 70:
                boot_x.append(pos[id(clade)]["x"])
                boot_y.append(pos[id(clade)]["y"])
                boot_text.append(str(val))
    if boot_x:
        fig.add_trace(go.Scatter(
            x=boot_x, y=boot_y,
            mode="text",
            text=boot_text,
            textposition="top center",
            textfont=dict(size=7, color="#333333"),
            hoverinfo="skip",
            showlegend=False,
            xaxis="x1", yaxis="y1",
        ))

    # ── 2. Leaf labels (own narrow column so they don't clip into the tree) ──
    fig.add_trace(go.Scatter(
        x=[0.0] * n, y=ly,
        mode="text",
        text=names,
        textposition="middle right",
        textfont=dict(size=8, color="#222222"),
        hovertext=act_hover,
        hoverinfo="text",
        showlegend=False,
        xaxis="x5", yaxis="y1",   # dedicated label axis
    ))

    # ── 3. Actuation bars ─────────────────────────────────────────────────────
    fig.add_trace(go.Bar(
        x=[1] * n,
        y=ly,
        orientation="h",
        marker=dict(color=colors, line=dict(color="#aaaaaa", width=0.4)),
        width=0.82,
        hovertext=act_hover,
        hoverinfo="text",
        showlegend=False,
        xaxis="x2", yaxis="y1",
    ))

    # ── 4. Legend dummy traces for actuation bins ─────────────────────────────
    for color, label in zip(ACT_COLORS, ACT_LABELS):
        fig.add_trace(go.Bar(
            x=[None], y=[None],
            orientation="h",
            marker=dict(color=color, line=dict(color="#aaaaaa", width=0.4)),
            name=label,
            legendgroup="act",
            showlegend=True,
        ))
    fig.add_trace(go.Bar(
        x=[None], y=[None],
        orientation="h",
        marker=dict(color=ACT_MISSING, line=dict(color="#aaaaaa", width=0.4)),
        name="missing",
        legendgroup="act",
        showlegend=True,
    ))

    # ── 5. MSA heatmap ────────────────────────────────────────────────────────
    msa_trace_start = len(fig.data)
    fig.add_trace(go.Heatmap(
        z=msa_z,
        y=ly,
        x=aln_positions,
        text=msa_cell,
        texttemplate="%{text}",
        textfont=dict(size=6, color="black"),
        customdata=msa_hover,
        hovertemplate="%{customdata}<extra></extra>",
        colorscale=DNA_COLORSCALE,
        zmin=0, zmax=5,
        zauto=False,
        showscale=False,
        showlegend=False,
        ygap=1,
        xaxis="x3", yaxis="y1",
    ))

    # ── 6. Sequence logo (stacked frequency bars) ─────────────────────────────
    logo_base = [0.0] * aln_len
    for nuc in "ATGC":
        f = freqs[nuc]
        fig.add_trace(go.Bar(
            x=aln_positions,
            y=f,
            base=logo_base,
            marker=dict(color=DNA_COLORS[nuc], line=dict(width=0)),
            name=nuc,
            legendgroup="dna",
            showlegend=False,
            hovertemplate=f"pos %{{x}}: {nuc} {nuc} %{{y:.1%}}<extra></extra>",
            xaxis="x4", yaxis="y2",
        ))
        logo_base = [a + b for a, b in zip(logo_base, f)]

    msa_trace_end = len(fig.data)
    msa_indices   = list(range(msa_trace_start, msa_trace_end))
    n_traces      = msa_trace_end
    vis_shown  = [True]  * n_traces
    vis_hidden = [False if i in set(msa_indices) else True for i in range(n_traces)]

    # ── Layout ────────────────────────────────────────────────────────────────
    # Row height adapts to keep the figure landscape (wider than tall).
    # Base width is 1400px; row height grows for small n, is floored at 8px.
    fig_width = 1400
    row_h     = max(8, min(16, int(fig_width * 0.58 / max(n, 1))))
    fig_height = max(350, n * row_h + 110)
    # If still portrait, widen until landscape
    if fig_height > fig_width:
        fig_width = int(fig_height * 1.45)

    label_fs = max(8, min(11, row_h))   # leaf label font scales with row height

    fig.update_layout(
        title=dict(text=peak_id, font=dict(size=13, family="Arial"), x=0.0, xanchor="left"),
        height=fig_height,
        width=fig_width,
        plot_bgcolor="white",
        paper_bgcolor="white",
        bargap=0,
        barmode="overlay",
        # Wipe out the default Plotly template so its heatmap colorscale
        # does not override the explicit colorscale on the MSA trace.
        template="none",
        margin=dict(l=10, r=140, t=45, b=40),
        font=dict(family="Arial"),
        legend=dict(
            title=dict(text="FIRE actuation", font=dict(size=10)),
            x=1.01, y=0.5,
            xanchor="left", yanchor="middle",
            tracegroupgap=0,
            font=dict(size=10),
        ),

        # ── Axes ──────────────────────────────────────────────────────────────
        xaxis1=dict(
            domain=TREE_X,
            anchor="y1",
            range=[-0.01 * max_depth, max_depth * 1.05],
            title=dict(text="branch length", font=dict(size=10)),
            showgrid=False, zeroline=True, zerolinecolor="#dddddd",
            tickfont=dict(size=9),
        ),
        xaxis5=dict(
            domain=LABEL_X,
            anchor="y1",
            range=[0, 1],
            showticklabels=False, showgrid=False, zeroline=False,
        ),
        xaxis2=dict(
            domain=ACT_X,
            anchor="y1",
            title=dict(text="% FIRE", font=dict(size=10)),
            showticklabels=False, showgrid=False, zeroline=False,
            range=[0, 1.1],
        ),
        xaxis3=dict(
            domain=MSA_X,
            anchor="y1",
            title=dict(text="alignment position", font=dict(size=10)),
            showgrid=False, zeroline=False,
            tickfont=dict(size=9),
        ),
        xaxis4=dict(
            domain=MSA_X,
            anchor="y2",
            showticklabels=False,
            showgrid=False, zeroline=False,
            matches="x3",   # zoom/pan stays in sync with MSA x-axis
        ),
        yaxis1=dict(
            domain=MAIN_Y,
            range=Y_RANGE,
            showticklabels=False, showgrid=False, zeroline=False,
        ),
        yaxis2=dict(
            domain=LOGO_Y,
            anchor="x4",
            range=[0, 1.05],
            title=dict(text="freq.", font=dict(size=9)),
            showgrid=False, zeroline=False, tickfont=dict(size=8),
        ),

        # ── MSA toggle: uses method="update" so domains are reclaimed ─────────
        updatemenus=[dict(
            type="buttons",
            direction="right",
            x=MSA_X[0],
            y=1.04,
            xanchor="left",
            buttons=[
                dict(
                    label="Show MSA",
                    method="update",
                    args=[
                        {"visible": vis_shown},
                        {
                            "xaxis1.domain": TREE_X,
                            "xaxis5.domain": LABEL_X,
                            "xaxis2.domain": ACT_X,
                            "xaxis3.visible": True,
                            "xaxis4.visible": True,
                            "yaxis2.visible": True,
                        },
                    ],
                ),
                dict(
                    label="Hide MSA",
                    method="update",
                    args=[
                        {"visible": vis_hidden},
                        {
                            "xaxis1.domain": TREE_X_H,
                            "xaxis5.domain": LABEL_X_H,
                            "xaxis2.domain": ACT_X_H,
                            "xaxis3.visible": False,
                            "xaxis4.visible": False,
                            "yaxis2.visible": False,
                        },
                    ],
                ),
            ],
            active=0,   # "Show MSA" is default
            showactive=True,
            bgcolor="#f4f4f4",
            bordercolor="#cccccc",
            font=dict(size=10),
        )],
    )

    # Update leaf label font to the adaptive size (trace was added earlier)
    for trace in fig.data:
        if getattr(trace, "xaxis", None) == "x5":
            trace.textfont = dict(size=label_fs, color="#222222")

    # DNA nucleotide legend — placed to the right of the logo panel
    for i, nuc in enumerate("ATGC"):
        fig.add_annotation(
            text=f"<b>{nuc}</b>",
            x=1.01,
            y=LOGO_Y[1] - (i + 0.5) * (LOGO_Y[1] - LOGO_Y[0]) / 4,
            xref="paper", yref="paper",
            showarrow=False,
            font=dict(size=10, color=DNA_COLORS[nuc]),
            xanchor="left",
        )

    return fig


# ── CLI ───────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Interactive tree + FIRE actuation + MSA plot.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--tree",      required=True, type=Path)
    parser.add_argument("--actuation", required=True, type=Path)
    parser.add_argument("--alignment", required=True, type=Path,
                        help="MAFFT-aligned FASTA (aligned/{peak_id}.fa)")
    parser.add_argument("--html",      required=True, type=Path,
                        help="Output HTML path")
    parser.add_argument("--pdf",       required=True, type=Path,
                        help="Output PDF path (via kaleido)")
    parser.add_argument("--peak-id",   type=str, default=None)
    args = parser.parse_args()

    tree               = Phylo.read(str(args.tree), "newick")
    actuation, coverage = load_actuation(args.actuation)
    alignment          = load_alignment(args.alignment)
    peak_id            = args.peak_id or args.tree.stem

    fig = make_figure(tree, actuation, coverage, alignment, peak_id)

    args.html.parent.mkdir(parents=True, exist_ok=True)

    # HTML – self-contained, toolbar "Download" button saves SVG
    fig.write_html(
        str(args.html),
        include_plotlyjs="cdn",
        config={
            "toImageButtonOptions": {"format": "svg", "filename": peak_id},
            "displayModeBar": True,
        },
    )
    print(f"Saved {args.html}", file=sys.stderr)

    # PDF via kaleido
    try:
        fig.write_image(str(args.pdf), format="pdf", scale=2)
        print(f"Saved {args.pdf}", file=sys.stderr)
    except Exception as e:
        print(f"WARNING: PDF export failed ({e}). Is kaleido installed?", file=sys.stderr)


if __name__ == "__main__":
    main()
