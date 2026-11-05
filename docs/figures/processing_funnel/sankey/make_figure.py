"""Processing funnel as a left-to-right Sankey (matplotlib Bezier ribbons).

HSCA peaks -> haplotype-selective -> heritable -> peaks with a driver ->
peaks with drivers on GRCh38 (SNV count omitted: inflated by linked blocks).

Computes counts from the run tables (../funnel_counts.py) and writes funnel.{svg,pdf,png} next to this script.
Run from the repo root:

    pixi run python docs/figures/processing_funnel/sankey/make_figure.py

Geometry is in millimetres. Node bar height is proportional to sqrt(count) so
the 3,234-peak end of the funnel stays readable next to 180,684 HSCA peaks.
At each node the continuing flow takes the next node's height and the dropped
categories share the remainder in proportion to their counts.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import PathPatch, Rectangle  # noqa: E402
from matplotlib.path import Path as MPath  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from funnel_counts import counts  # noqa: E402

C = counts(sys.argv[1] if len(sys.argv) > 1 else "hsca_hapsel")
OUT = HERE / "funnel"

# ------------------------------------------------------------------ style
INK, GREY, MIDGREY = "#222222", "#6b6b6b", "#9a9a9a"
TRUNK, DROP = "#4a4a4a", "#dcdcdc"
RED, DARKRED = "#d62728", "#a50f15"
RED_RIB = RED
plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

W_MM, H_MM = 120.0, 64.0
K = 36.0 / C["hsca"] ** 0.5  # mm per sqrt(count): HSCA bar is 36 mm tall
NODE_W = 1.3  # node bar width (mm)
MIN_DROP = 0.45  # thinnest drawn ribbon (mm)
GAP = 1.4  # vertical gap between lanes (mm)
LABEL_H = 2.6  # vertical room a one-line 6 pt label needs (mm)
Y_TOP = 49.5  # flat top edge of the main path


def h(n: int) -> float:
    """Bar height in mm for a count (sqrt scale)."""
    return K * n**0.5


def fmt(n: int) -> str:
    return f"{n:,}"


# ------------------------------------------------------------------ geometry
def ribbon(ax, x0, t0, b0, x1, t1, b1, color, z=1, lw=0.0, ec=None):
    """Ribbon from x0 [t0..b0] to x1 [t1..b1]; edges are cubic Beziers."""
    xm = (x0 + x1) / 2
    verts = [
        (x0, t0), (xm, t0), (xm, t1), (x1, t1),
        (x1, b1),
        (xm, b1), (xm, b0), (x0, b0),
        (x0, t0),
    ]
    codes = [MPath.MOVETO] + [MPath.CURVE4] * 3 + [MPath.LINETO] \
        + [MPath.CURVE4] * 3 + [MPath.CLOSEPOLY]
    ax.add_patch(PathPatch(MPath(verts, codes), fc=color, ec=ec or color, lw=lw,
                           zorder=z, joinstyle="round"))


def pair_label(ax, x, y, bold: str, rest: str, size: float = 6) -> float:
    """Bold count followed by a regular label; returns the right edge (mm)."""
    t = ax.text(x, y, bold, ha="left", va="center", fontsize=size,
                color=GREY, fontweight="bold")
    r = ax.figure.canvas.get_renderer()
    bb = t.get_window_extent(r).transformed(ax.transData.inverted())
    t2 = ax.text(bb.x1 + 0.9, y, rest, ha="left", va="center",
                 fontsize=size, color=GREY)
    return t2.get_window_extent(r).transformed(ax.transData.inverted()).x1


@dataclass
class Node:
    x: float
    n: int
    count_label: str
    lines: list[str]
    color: str = INK

    @property
    def h(self) -> float:
        return h(self.n)


@dataclass
class Drop:
    src: int  # index of source node
    n: int
    text: str
    level: int  # lane level below the main path (0 = nearest)
    run: float  # horizontal run before the label (mm)
    width: float = 0.0
    y_src_top: float = 0.0


def main() -> None:
    nodes = [
        Node(8.0, C["hsca"], fmt(C["hsca"]), ["HSCA peaks"]),
        Node(32.0, C["hapsel"], fmt(C["hapsel"]),
             ["haplotype-selective", "peaks (analyzed)"]),
        Node(56.5, C["heritable"], fmt(C["heritable"]), ["heritable", "peaks"]),
        Node(81.0, C["drv_peaks"], fmt(C["drv_peaks"]),
             ["peaks with a driver", f"{fmt(C['drv_sets'])} driver sets"], RED),
        Node(107.5, C["drv_peaks_hg38"], fmt(C["drv_peaks_hg38"]),
             ["peaks with drivers", "on GRCh38"], RED),
    ]
    # Drops from the same node travel as one bundle into a shared lane.
    drops = [
        Drop(0, C["epi"], "epigenetically selective", 2, 11.0),
        Drop(0, C["gen"], "genetically specific (indels)", 2, 11.0),
        Drop(0, C["complex"], "complex", 2, 11.0),
        Drop(1, C["not_heritable"], "not heritable", 1, 9.0),
        Drop(2, C["her_no_drv"], "no driver", 0, 6.0),
        # peaks whose only drivers are indels have no SNV to place on GRCh38
        Drop(3, C["drv_peaks"] - C["drv_peaks_hg38"], "not placed on GRCh38", 0, 6.0),
    ]

    # Split each node: continuing flow on top, drops underneath.
    keep_src = {}  # node index -> height of continuing flow at the source
    for i, nd in enumerate(nodes[:-1]):
        ds = [d for d in drops if d.src == i]
        nxt = nodes[i + 1].h
        total = sum(d.n for d in ds)
        room = max(nd.h - nxt, MIN_DROP * len(ds)) if ds else 0.0
        for d in ds:
            d.width = max(room * d.n / total, MIN_DROP)
        room = sum(d.width for d in ds)
        keep_src[i] = nd.h - room
        y = Y_TOP - keep_src[i]
        for d in ds:  # stacked below the continuing flow, in list order
            d.y_src_top = y
            y -= d.width

    fig = plt.figure(figsize=(W_MM / 25.4, H_MM / 25.4))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W_MM)
    ax.set_ylim(0, H_MM)
    ax.set_aspect("equal")
    ax.axis("off")

    # ---- main path ribbons
    for i in range(len(nodes) - 1):
        a, b = nodes[i], nodes[i + 1]
        col = RED_RIB if b.color == RED else TRUNK
        ribbon(ax, a.x + NODE_W / 2, Y_TOP, Y_TOP - keep_src[i],
               b.x - NODE_W / 2, Y_TOP, Y_TOP - b.h, col, z=2, lw=0.2)

    # ---- drop lanes: one lane per level, later stages nearest the path
    def bundle_height(lvl: int) -> float:
        ds = [d for d in drops if d.level == lvl]
        span = sum(d.width for d in ds)
        return max(span, LABEL_H * len(ds) - (ds[-1].width if len(ds) > 1 else 0))

    lane_top = {}
    y = Y_TOP - max(nd.h for nd in nodes[2:]) - GAP
    for lvl in sorted({d.level for d in drops}):
        lane_top[lvl] = y
        y -= bundle_height(lvl) + GAP
    last_label_y: dict[int, float] = {}
    label_y: dict[str, float] = {}
    note_x = 60.0
    for d in drops:
        x0 = nodes[d.src].x + NODE_W / 2
        group = [g for g in drops if g.src == d.src]
        span = sum(g.width for g in group)
        drop_by = group[0].y_src_top - lane_top[d.level]
        if len(group) == 1 and d.width < LABEL_H:
            drop_by += (LABEL_H - d.width) / 2  # centre thin ribbon in its lane
        drop_by = max(drop_by, 0.0)
        t1 = d.y_src_top - drop_by
        x_bend = x0 + d.run
        ribbon(ax, x0, d.y_src_top, d.y_src_top - d.width,
               x_bend, t1, t1 - d.width, DROP, z=1, lw=0.35, ec="white")
        yc = t1 - d.width / 2
        if d.src in last_label_y:  # keep stacked labels apart
            yc = min(yc, last_label_y[d.src] - LABEL_H)
        last_label_y[d.src] = yc
        label_y[d.text] = yc
        if d.src == 2:
            note_x = x_bend + 1.0
        pair_label(ax, x_bend + 1.0, yc, fmt(d.n), d.text)
        del span

    # ---- node bars and labels
    for nd in nodes:
        dark = DARKRED if nd.color == RED else INK
        ax.add_patch(Rectangle((nd.x - NODE_W / 2, Y_TOP - nd.h), NODE_W, nd.h,
                               fc=dark, ec="none", zorder=3))
        ax.text(nd.x, Y_TOP + 9.6, nd.count_label, ha="center", va="baseline",
                fontsize=9, fontweight="bold", color=INK)
        for j, line in enumerate(nd.lines):
            ax.text(nd.x, Y_TOP + 6.6 - 2.3 * j, line, ha="center",
                    va="baseline", fontsize=6, color=GREY)

    ax.text(note_x, label_y["complex"], "bar height scales with \u221acount",
            ha="left", va="center", fontsize=5.5, color=MIDGREY)

    for ext in ("svg", "pdf"):
        fig.savefig(f"{OUT}.{ext}")
    fig.savefig(f"{OUT}.png", dpi=300)


if __name__ == "__main__":
    main()
