"""Schematic of per-FIRE-peak haplotype classes (panel e).

Run from the repo root:
    pixi run python docs/figures/hsca_schematic/make_figure.py

Writes hsca_schematic.{svg,pdf,png} next to this script. All coordinates
are in millimetres on a 120 x 46 mm canvas.
"""

from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Rectangle  # noqa: E402

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

OUT = Path(__file__).resolve().parent
W_MM, H_MM = 120.0, 46.0

# colours
PAT_LINE, PAT_BAR, PAT_BLOCK, PAT_TEXT = "#fa6800", "#f8c9a8", "#f2a46e", "#c75300"
MAT_LINE, MAT_BAR, MAT_BLOCK, MAT_TEXT = "#5b9bd5", "#bcd3ee", "#8fb5e3", "#2f6fae"
FIRE = "#8b0000"
INK = "#222222"
GUIDE = "#7f7f7f"

FS = 6.0  # base font size (pt)
FS_SMALL = 5.0

# layout (mm)
LEFT = 13.0  # label column
COL_W = (W_MM - LEFT - 0.5) / 5
HALF = 7.4  # half-width of a locus
N_FIB = 6
FIB_DY = 1.45
TICK_W, TICK_H = 1.6, 0.85
BAR_H = 1.9
PEAK_MS = 3.4  # FIRE peak triangle marker size (pt)

Y_HEADER = 44.7  # HSCA family label
Y_BRACKET = 43.1
PAT_FIB_Y = [41.0 - i * FIB_DY for i in range(N_FIB)]
PAT_BAR_Y = 30.9
BOX_Y, BOX_H, BOX_W = 26.4, 3.6, 20.0
DESC_Y = 22.2
MAT_BAR_Y = 17.9
MAT_FIB_Y = [14.8 - i * FIB_DY for i in range(N_FIB)]
Y_FOOTER = 2.2  # FIRE call-out and FIRE peak key

# per-fiber start/end jitter so fibers look like reads, not a grid
JITTER = [(-0.4, 0.6), (0.9, -0.2), (-0.1, 1.1), (0.5, 0.3), (1.2, -0.5), (0.0, 0.8)]

HIGH, LOW = [0, 1, 2, 3, 5], [2]  # 5/6 vs 1/6 fibers with a FIRE

CLASSES = [  # one entry per drawn locus (HSCA-SV has two)
    {
        "pat_re": [0.0],
        "mat_re": [0.0],
        "links": [(0.0, 0.0)],
        "pat_fire": [[0, 1, 2, 4, 5]],
        "mat_fire": [[0, 1, 3, 4, 5]],
        "frac": ("5/6", "5/6"),
    },
    {
        "pat_re": [0.0],
        "mat_re": [0.0],
        "links": [(0.0, 0.0)],
        "pat_fire": [HIGH],
        "mat_fire": [LOW],
        "frac": ("5/6", "1/6"),
        "bases": ("A", "A", "#555555"),
    },
    {
        "pat_re": [0.0],
        "mat_re": [0.0],
        "links": [(0.0, 0.0)],
        "pat_fire": [HIGH],
        "mat_fire": [LOW],
        "frac": ("5/6", "1/6"),
        "bases": ("A", "G", "#8b0000"),
    },
    {
        "pat_re": [0.0],
        "mat_re": [],
        "links": [(0.0, None)],
        "pat_fire": [[0, 1, 2, 3, 4]],
        "mat_fire": [],
        "pat_blocks": [(0.0, 5.6)],  # (centre offset, width): the insertion
        "mat_gap": 2.8,  # half-width of the matching gap on the maternal side
    },
    {
        "pat_re": [-3.7, 0.0, 3.7],
        "mat_re": [0.0],
        "links": [(-3.7, 0.0), (0.0, 0.0), (3.7, 0.0)],
        "pat_fire": [[0, 1, 3, 4], [1, 2, 4, 5], [0, 2, 3, 5]],
        "mat_fire": [[0, 1, 2, 3, 5]],
        "pat_blocks": [(-3.7, 3.3), (0.0, 3.3), (3.7, 3.3)],  # tandem copies
        "mat_blocks": [(0.0, 3.3)],
    },
]

# class boxes: (name, descriptor, colour, loci covered). HSCA-SV spans two
# sub-examples: a FIRE peak on one haplotype only, and tandem copies collapsing.
BOXES = [
    ("Shared", "similar accessibility", "#888888", [0]),
    ("HSCA-noVar", "same sequence,\ndiffers epigenetically", "#800080", [1]),
    ("HSCA-var", "SNVs / small indels", "#8b0000", [2]),
    ("HSCA-SV", "no 1:1 homolog (SV)", "#f0a500", [3, 4]),
]


def segments(cx: float, gap: float | None) -> list[tuple[float, float]]:
    """Return the (x0, x1) spans of a locus, split around a central gap."""
    x0, x1 = cx - HALF, cx + HALF
    if gap is None:
        return [(x0, x1)]
    return [(x0, cx - gap), (cx + gap, x1)]


def draw_fibers(ax, cx, ys, color, re_pos, fire_sets, gap=None, shift=0):
    """Draw a stack of fibers with FIRE ticks over the FIRE peak positions."""
    for i, y in enumerate(ys):
        j0, j1 = JITTER[(i + shift) % len(JITTER)]
        segs = segments(cx, gap)
        segs[0] = (segs[0][0] + j0, segs[0][1])
        segs[-1] = (segs[-1][0], segs[-1][1] - j1)
        for a, b in segs:
            ax.plot([a, b], [y, y], color=color, lw=0.6, solid_capstyle="butt", zorder=2)
        for pos, fired in zip(re_pos, fire_sets):
            if i in fired:
                ax.add_patch(
                    Rectangle(
                        (cx + pos - TICK_W / 2, y - TICK_H / 2),
                        TICK_W,
                        TICK_H,
                        facecolor=FIRE,
                        edgecolor="none",
                        zorder=3,
                    )
                )


def draw_bar(ax, cx, y, color, block_color, re_pos, blocks=(), gap=None):
    """Draw an assembly segment, optional peak blocks, a gap, and FIRE peak triangles."""
    for a, b in segments(cx, gap):
        ax.add_patch(
            Rectangle(
                (a, y - BAR_H / 2), b - a, BAR_H, facecolor=color, edgecolor="none", zorder=2
            )
        )
    if gap is not None:
        # alignment gap: the haplotype skips this sequence
        ax.plot(
            [cx - gap, cx + gap], [y, y], color=GUIDE, lw=0.5, ls=(0, (1.0, 1.0)), zorder=2
        )
    for off, w in blocks:
        ax.add_patch(
            Rectangle(
                (cx + off - w / 2, y - BAR_H / 2),
                w,
                BAR_H,
                facecolor=block_color,
                edgecolor="white",
                lw=0.4,
                zorder=3,
            )
        )
    for pos in re_pos:
        ax.plot(cx + pos, y, marker="^", ms=PEAK_MS, mfc=FIRE, mec="none", zorder=4)


def main() -> None:
    fig = plt.figure(figsize=(W_MM / 25.4, H_MM / 25.4))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W_MM)
    ax.set_ylim(0, H_MM)
    ax.axis("off")

    centers = [LEFT + COL_W * (i + 0.5) for i in range(len(CLASSES))]

    for k, (cx, c) in enumerate(zip(centers, CLASSES)):
        gap = c.get("mat_gap")
        draw_fibers(ax, cx, PAT_FIB_Y, PAT_LINE, c["pat_re"], c["pat_fire"], shift=k)
        draw_fibers(
            ax, cx, MAT_FIB_Y, MAT_LINE, c["mat_re"], c["mat_fire"], gap=gap, shift=k + 3
        )
        draw_bar(ax, cx, PAT_BAR_Y, PAT_BAR, PAT_BLOCK, c["pat_re"], c.get("pat_blocks", ()))
        draw_bar(
            ax, cx, MAT_BAR_Y, MAT_BAR, MAT_BLOCK, c["mat_re"], c.get("mat_blocks", ()), gap
        )

        # dashed homology links between the two assemblies
        y_top = PAT_BAR_Y - BAR_H / 2
        y_bot = MAT_BAR_Y + (BAR_H / 2 if gap is None else 0.0)
        for p, m in c["links"]:
            xm = cx + (0.0 if m is None else m)
            ax.plot(
                [cx + p, xm],
                [y_top, y_bot],
                color=GUIDE,
                lw=0.5,
                ls=(0, (2.2, 1.6)),
                zorder=1,
            )

        # sequence at the FIRE peak: identical (noVar) or a SNV (var)
        if "bases" in c:
            pa, ma, col = c["bases"]
            for y, base in ((PAT_BAR_Y, pa), (MAT_BAR_Y, ma)):
                ax.text(
                    cx + 3.0,
                    y,
                    base,
                    ha="center",
                    va="center_baseline",
                    fontsize=5.2,
                    fontweight="bold",
                    color=col,
                    zorder=5,
                )

        # FIRE fraction next to each fiber stack
        if "frac" in c:
            for ys, txt in ((PAT_FIB_Y, c["frac"][0]), (MAT_FIB_Y, c["frac"][1])):
                ax.text(
                    cx + HALF + 0.9,
                    sum(ys) / N_FIB,
                    txt,
                    ha="left",
                    va="center",
                    fontsize=FS_SMALL,
                    color=FIRE,
                    zorder=5,
                )

    for name, desc, color, loci in BOXES:
        cx = sum(centers[i] for i in loci) / len(loci)
        bw = BOX_W + COL_W * (len(loci) - 1)
        ax.add_patch(
            FancyBboxPatch(
                (cx - bw / 2, BOX_Y - BOX_H / 2),
                bw,
                BOX_H,
                boxstyle="round,pad=0,rounding_size=0.6",
                facecolor=color,
                edgecolor="none",
                zorder=5,
            )
        )
        ax.text(
            cx,
            BOX_Y,
            name,
            ha="center",
            va="center",
            fontsize=FS,
            fontweight="bold",
            color="white",
            zorder=6,
        )
        ax.text(
            cx,
            DESC_Y,
            desc,
            ha="center",
            va="center",
            fontsize=FS_SMALL,
            color="#444444",
            linespacing=1.05,
            zorder=6,
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.6},
        )

    # HSCA family bracket over columns 2-5
    bx0, bx1 = centers[1] - BOX_W / 2, centers[-1] + BOX_W / 2
    ax.plot(
        [bx0, bx0, bx1, bx1],
        [Y_BRACKET - 0.9, Y_BRACKET, Y_BRACKET, Y_BRACKET - 0.9],
        color=INK,
        lw=0.5,
        solid_joinstyle="miter",
    )
    ax.text(
        (bx0 + bx1) / 2,
        Y_HEADER,
        "HSCA: haplotype-selective chromatin accessibility",
        ha="center",
        va="center",
        fontsize=FS,
        fontweight="bold",
        color=INK,
    )

    # row labels
    def row_label(y, text, color):
        ax.text(0.6, y, text, ha="left", va="center", fontsize=FS, color=color, linespacing=1.0)

    row_label(sum(PAT_FIB_Y) / N_FIB, "Paternal\nfibers", PAT_TEXT)
    row_label(PAT_BAR_Y, "Paternal\nassembly", PAT_TEXT)
    row_label(MAT_BAR_Y, "Maternal\nassembly", MAT_TEXT)
    row_label(sum(MAT_FIB_Y) / N_FIB, "Maternal\nfibers", MAT_TEXT)

    # FIRE call-out with a dotted leader from the bottom tick in the first column
    tick_x, tick_y = centers[0], MAT_FIB_Y[-1] - TICK_H / 2
    lx = centers[0] + 6.0
    ax.plot(
        [tick_x, tick_x, lx - 0.6],
        [tick_y - 0.3, Y_FOOTER, Y_FOOTER],
        color=INK,
        lw=0.45,
        ls=(0, (0.6, 0.9)),
        zorder=4,
    )
    callout = ax.text(
        lx,
        Y_FOOTER,
        "Fiber-seq inferred regulatory element (FIRE)",
        ha="left",
        va="center",
        fontsize=FS,
        color=FIRE,
    )

    # FIRE peak key, ~6 mm after the call-out
    r = fig.canvas.get_renderer()
    sx = callout.get_window_extent(r).transformed(ax.transData.inverted()).x1 + 6.0
    ax.plot(sx, Y_FOOTER, marker="^", ms=PEAK_MS, mfc=FIRE, mec="none")
    ax.text(sx + 1.4, Y_FOOTER, "FIRE peaks", ha="left", va="center", fontsize=FS, color=INK)

    for ext, kw in (("svg", {}), ("pdf", {}), ("png", {"dpi": 300})):
        fig.savefig(OUT / f"hsca_schematic.{ext}", facecolor="white", **kw)
    plt.close(fig)


if __name__ == "__main__":
    main()
