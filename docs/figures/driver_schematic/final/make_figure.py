#!/usr/bin/env python
"""Driver-variant schematic panel, final ("one tree, built up").

Built from the exemplar peak chr8_70445067_70445271 (../toy_data.json): the
14-tip subsampled tree, per-tip FIRE actuation, carriers of RE112 (closing)
and RE63 (opening), the H2 estimate, the forward-selection R2 steps and the
joint-model betas/P values. The same tree is reused across steps so the
reader watches one object gain information.

Illustrative only (not in toy_data.json): the horizontal offsets and the
small flank indels that make step 1 look unaligned, and the matching gap
marks in step 2. The real exemplar alignment is gap-free (all 78 sequences
704 bp), so these show what alignment does in general.

Run from the repo root:
    pixi run python docs/figures/driver_schematic/final/make_figure.py
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Polygon, Rectangle  # noqa: E402

N_INDIVIDUALS = 39

HERE = Path(__file__).resolve().parent
OUT = HERE / "driver_schematic"
DATA = json.loads((HERE.parent / "toy_data.json").read_text())

# ------------------------------------------------------------------ style
W = 120.0  # mm; height is set from the layout below
MM = 1 / 25.4
FS_TITLE, FS_SUB, FS_SMALL = 5.6, 6.0, 5.5
LW = 0.6
# 10 actuation bins (0-10% ... 90-100%), as in plot_tree_compact.py
ACT_COLORS = ["#ffffff", "#fff5f0", "#fee0d2", "#fcbba1", "#fc9272",
              "#fb6a4a", "#ef3b2c", "#cb181d", "#a50f15", "#67000d"]
CLOSE, OPEN = "#d62728", "#1f77b4"
INK, DARK, GREY, LIGHT, FAINT = "#222222", "#444444", "#6b6b6b", "#bbbbbb", "#e4e4e4"
VAR_GREY, GAP = "#8a8a8a", "#555555"

TITLES = [
    f"Extract FIRE peaks ± 250 bp from {2 * N_INDIVIDUALS} haplotypes",
    "Align (MAFFT), build tree (IQ-TREE)",
    "Keep heritable peaks (MCMCglmm)",
    "Forward-select variants (ΔR²\u00a0≥\u00a05%)",
    "Drivers: joint β and P",
]
# title boxes (x0, x1) in mm; icons below are centred on these
TITLE_BOX = [(0.0, 23.5), (23.5, 42.5), (42.5, 66.0), (66.0, 95.0), (95.0, 120.0)]
LABEL = {"RE112": "Driver 1", "RE63": "Driver 2"}  # display names, in selection order
H2_DEF = "H²ₚ: phylogenetic heritability"
R2_DEF = "R²: variance explained by variants (coverage-weighted regression)"
FS_DEF = 5.0

plt.rcParams.update({
    "font.family": ["Arial", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "lines.solid_capstyle": "butt",
})


# ------------------------------------------------------------------ title wrapping
def wrap(s: str, width_mm: float, size: float, weight: str = "bold") -> list[str]:
    """Greedy word wrap using real rendered text widths."""
    f = plt.figure(figsize=(4, 1))
    r = f.canvas.get_renderer()

    def wmm(t):
        tt = f.text(0, 0, t, fontsize=size, fontweight=weight)
        w = tt.get_window_extent(renderer=r).width / f.dpi * 25.4
        tt.remove()
        return w

    out, cur = [], ""
    for word in s.split(" "):  # split on plain spaces only; NBSP keeps phrases together
        trial = f"{cur} {word}".strip()
        if cur and wmm(trial) > width_mm:
            out.append(cur)
            cur = word
        else:
            cur = trial
    out.append(cur)
    plt.close(f)
    return out


TITLE_LINES = [wrap(t, x1 - x0 - 1.0, FS_TITLE) for t, (x0, x1) in zip(TITLES, TITLE_BOX)]
LINE_H = FS_TITLE * 1.22 / 72 * 25.4
TB = 0.8 + max(len(t) for t in TITLE_LINES) * LINE_H  # titles end here

# ------------------------------------------------------------------ vertical layout
TIPS = DATA["tips"]  # newick order, top to bottom
TIP_IDX = {t["tip"]: i for i, t in enumerate(TIPS)}
N = len(TIPS)
DY = 1.9
T0 = TB + 4.3
TIP_Y = [T0 + i * DY for i in range(N)]
YL = TIP_Y[-1]
MID = (TIP_Y[0] + YL) / 2
H = YL + 12.4

fig = plt.figure(figsize=(W * MM, H * MM))
ax = fig.add_axes((0, 0, 1, 1))
ax.set_xlim(0, W)
ax.set_ylim(H, 0)  # y grows downward, mm
ax.axis("off")


def act_color(a: float) -> str:
    """Bin colour for an actuation fraction in [0, 1]."""
    return ACT_COLORS[min(int(a * 10), 9)]


def line(xs, ys, color=INK, lw=LW, **kw):
    ax.plot(xs, ys, color=color, linewidth=lw, **kw)


def rect(x, y, w, h, fc, ec="none", lw=0.3, **kw):
    ax.add_patch(Rectangle((x, y), w, h, facecolor=fc, edgecolor=ec, linewidth=lw, **kw))


def text(x, y, s, size=FS_SUB, color=GREY, ha="center", va="center", **kw):
    ax.text(x, y, s, fontsize=size, color=color, ha=ha, va=va, **kw)


CELL_W, CELL_H = 1.5, 1.35


def cell(x, y, a):
    """Actuation cell, left edge at x, centred vertically on y."""
    rect(x, y - CELL_H / 2, CELL_W, CELL_H, act_color(a), ec=GREY, lw=0.3)


for lines_, (x0, x1) in zip(TITLE_LINES, TITLE_BOX):
    text((x0 + x1) / 2, TB, "\n".join(lines_), size=FS_TITLE, color=INK,
         fontweight="bold", va="bottom", linespacing=1.15)

# ------------------------------------------------------------------ tree
# Real topology from toy_data.json tree_newick, with its schematic branch
# lengths: node -> (parent, branch length, children). Tips are leaves.
NODES = {
    "root": (None, 0.0, ["ancestral", "RE63", "RE112"]),
    "ancestral": ("root", 0.3, ["HG02257_2", "HG002_1", "HG00673_1", "HG01106_1", "HG04160_1"]),
    "RE63": ("root", 1.0, ["HG02071_2", "HG02132_2", "HG00438_2", "HG01106_2"]),
    "RE112": ("root", 1.0, ["bg1", "bg2"]),
    "bg1": ("RE112", 0.3, ["HG01891_1", "HG02257_1"]),
    "bg2": ("RE112", 2.0, ["RE159", "HG03579_2"]),
    "RE159": ("bg2", 1.0, ["HG02622_1", "HG00735_1"]),
}
PARENT = {k: n for n, (_, _, kids) in NODES.items() for k in kids}
TREE_DEPTH = 4.7  # tips drawn aligned just past the deepest node (4.0)
assert set(TIP_IDX) == {k for k in PARENT if k not in NODES}


def depth(n):
    d = 0.0
    while n != "root":
        d += NODES[n][1] if n in NODES else 0.0
        n = PARENT[n]
    return d


def node_y(n):
    if n in TIP_IDX:
        return TIP_Y[TIP_IDX[n]]
    ys = [node_y(c) for c in NODES[n][2]]
    return (min(ys) + max(ys)) / 2


def clade_color(n):
    """Driver colour for branches inside the RE112 / RE63 clades, else None."""
    while n != "root":
        if n == "RE112":
            return CLOSE
        if n == "RE63":
            return OPEN
        n = PARENT[n]
    return None


def draw_tree(x0, width, drivers=False):
    """Rectangular phylogram, root at x0, tips aligned at x0 + width."""

    def nx(n):
        return x0 + width * (1.0 if n in TIP_IDX else depth(n) / TREE_DEPTH)

    line([x0 - 0.6, x0], [node_y("root")] * 2)  # root stub
    for name, (_, _, kids) in NODES.items():
        ys = [node_y(k) for k in kids]
        c = (clade_color(name) if drivers else None) or INK
        line([nx(name)] * 2, [min(ys), max(ys)], color=c)
        for k in kids:
            ck = (clade_color(k) if drivers else None) or INK
            w = 1.5 if drivers and k in ("RE112", "RE63") else LW
            line([nx(name), nx(k)], [node_y(k)] * 2, color=ck, lw=w)
    return x0 + width


def chevron(x, y=MID):
    ax.add_patch(Polygon([(x - 0.6, y - 1.1), (x + 0.6, y), (x - 0.6, y + 1.1)],
                         closed=False, fill=False, edgecolor=LIGHT, linewidth=0.8))


# ------------------------------------------------------------------ sequences
PEAK = DATA["peak"]
L_TOT, FL, RE_BP = PEAK["aligned_cols"], PEAK["flank_bp"], PEAK["re_bp"]
BAR_H = 0.62
# No indels: the real exemplar alignment is ungapped (all sequences 704 bp).
INDELS: dict[str, tuple[int, int]] = {}


def segments(tip):
    """Ungapped (start, end) column runs of a tip on the alignment."""
    if tip not in INDELS:
        return [(0, L_TOT)]
    s, n = INDELS[tip]
    return [(0, s), (s + n, L_TOT)]


def seq_bar(x, y, scale, tip, aligned):
    """Grey bar (pale flanks, darker RE). Aligned bars keep gap columns as a thin
    dark line; unaligned bars simply omit them (the sequence is shorter)."""
    if aligned:
        for s, e in segments(tip):
            rect(x + s * scale, y - BAR_H / 2, (e - s) * scale, BAR_H, FAINT)
        rect(x + FL * scale, y - BAR_H / 2, RE_BP * scale, BAR_H, LIGHT)
        if tip in INDELS:
            s, n = INDELS[tip]
            line([x + s * scale, x + (s + n) * scale], [y, y], color=GAP, lw=0.5)
        return x + FL * scale
    length = L_TOT - (INDELS[tip][1] if tip in INDELS else 0)
    rect(x, y - BAR_H / 2, length * scale, BAR_H, FAINT)
    gap_before_re = INDELS[tip][1] if tip in INDELS and INDELS[tip][0] < FL else 0
    re_x = x + (FL - gap_before_re) * scale
    rect(re_x, y - BAR_H / 2, RE_BP * scale, BAR_H, LIGHT)
    return re_x


TRACK_MAX = DY - BAR_H - 0.12


def track(re_x, y, scale, a):
    """bedGraph-like FIRE profile above a bar, peaking over the RE."""
    base = y - BAR_H / 2 - 0.05
    xs = np.linspace(re_x - 0.4 * RE_BP * scale, re_x + 1.4 * RE_BP * scale, 40)
    c = re_x + 0.5 * RE_BP * scale
    h = TRACK_MAX * min(a / 0.7, 1.0)
    ys = base - h * np.exp(-0.5 * ((xs - c) / (0.32 * RE_BP * scale)) ** 2)
    if h > 0.05:
        ax.fill_between(xs, base, ys, color=act_color(max(a, 0.2)), linewidth=0,
                        step="mid")
        ax.step(xs, ys, where="mid", color="#cb181d", linewidth=0.25)
    line([xs[0], xs[-1]], [base, base], color=GREY, lw=0.25)


# ------------------------------------------------------------------ 1. input
SC1 = 11.4 / L_TOT
order = [11, 2, 7, 0, 9, 5, 13, 3, 8, 1, 12, 6, 4, 10]  # shuffled, not tree order
offs = [0.0, 1.3, -0.6, 0.7, 1.6, -0.3, 0.9, -0.8, 0.4, 1.4, -0.5, 0.2, 1.1, -0.2]
BX = 2.0
for i, (t, o) in enumerate(zip(order, offs)):
    tip = TIPS[t]
    re_x = seq_bar(BX + o, TIP_Y[i], SC1, tip["tip"], aligned=False)
    track(re_x, TIP_Y[i], SC1, tip["actuation"])
    cell(BX + 14.6, TIP_Y[i], tip["actuation"])
# RE +/- 250 bp bracket over the first bar (a closed haplotype: flat track)
bx0 = BX + offs[0]
bre0, bre1 = bx0 + FL * SC1, bx0 + (FL + RE_BP) * SC1
bx1 = bx0 + L_TOT * SC1
by = TIP_Y[0] - 1.25
line([bx0, bx0, bx1, bx1], [by + 0.4, by, by, by + 0.4], color=GREY, lw=0.4)
line([bre0, bre0], [by, by + 0.4], color=GREY, lw=0.4)
line([bre1, bre1], [by, by + 0.4], color=GREY, lw=0.4)
text((bx0 + bx1) / 2, by - 1.0, f"RE ± {FL} bp", size=FS_SMALL)
# key: FIRE actuation = track height/colour and square, 0 to 100%
ky = YL + 2.6
kx, kw = 4.2, 0.95
text(kx - 0.5, ky, "0", size=FS_SMALL, ha="right")
for b, c in enumerate(ACT_COLORS):
    rect(kx + b * kw, ky - 0.6, kw, 1.2, c)
rect(kx, ky - 0.6, 10 * kw, 1.2, "none", ec=GREY, lw=0.3)
text(kx + 10 * kw + 0.5, ky, "100%", size=FS_SMALL, ha="left")
text(kx + 5 * kw, ky + 2.2, "FIRE actuation", size=FS_SMALL)

# ------------------------------------------------------------------ 2. align + tree
tx = draw_tree(23.8, 5.6)
ABX, SC2 = tx + 0.8, 12.6 / L_TOT
# real alignment columns (0-based) of RE variable sites in the exemplar MSA
SITE_COL = {"RE63": 312, "RE74": 323, "RE112": 361, "RE148": 397, "RE159": 408}
VAR_SITES = [("RE74", "RE74_RE148", VAR_GREY), ("RE148", "RE74_RE148", VAR_GREY),
             ("RE159", "RE159", VAR_GREY), ("RE63", "RE63", OPEN),
             ("RE112", "RE112", CLOSE)]


def col_x(site):
    return ABX + (SITE_COL[site] + 0.5) * SC2


for i, t in enumerate(TIPS):
    seq_bar(ABX, TIP_Y[i], SC2, t["tip"], aligned=True)
    for site, key, c in VAR_SITES:
        if t[key]:
            line([col_x(site)] * 2, [TIP_Y[i] - 0.5, TIP_Y[i] + 0.5], color=c,
                 lw=0.7 if c != VAR_GREY else 0.5)
ly = YL + 1.9
text(col_x("RE63") - 0.25, ly, LABEL["RE63"], size=FS_SMALL, color=OPEN, ha="right")
text(col_x("RE112") + 0.25, ly, LABEL["RE112"], size=FS_SMALL, color=CLOSE, ha="left")
for s in ("RE63", "RE112"):
    c = OPEN if s == "RE63" else CLOSE
    line([col_x(s)] * 2, [YL + 0.6, ly - 0.1], color=c, lw=0.3)

# ------------------------------------------------------------------ 3. heritability
tx = draw_tree(47.8, 8.8)
for i, t in enumerate(TIPS):
    cell(tx + 0.6, TIP_Y[i], t["actuation"])
h2 = DATA["h2_examples"]["pass"]
hx = 54.4
text(hx, YL + 2.4, f"H²ₚ = {h2['h2']:.2f} [{h2['hpd'][0]:.2f}, {h2['hpd'][1]:.2f}]",
     color=DARK, size=FS_SMALL)
text(hx, YL + 3.6, "\n".join(wrap(H2_DEF, 20.5, FS_DEF, "normal")), color=GREY,
     size=FS_DEF, va="top", linespacing=1.15)
# gate: non-heritable peaks stop here
# ponytail: hand-placed stop sign; y positions tuned to clear the H² and R² notes
gx, gy = 64.4, YL + 8.6
ax.annotate("", xy=(gx, gy - 1.3), xytext=(gx, MID + 2.0),
            arrowprops=dict(arrowstyle="-|>", color=DARK, lw=0.8,
                            mutation_scale=5, shrinkA=0, shrinkB=0))
ax.add_patch(plt.Circle((gx, gy), 1.2, facecolor="white", edgecolor=DARK, lw=0.8, zorder=5))
line([gx - 0.6, gx + 0.6], [gy - 0.6, gy + 0.6], color=DARK, lw=0.8, zorder=6)
line([gx - 0.6, gx + 0.6], [gy + 0.6, gy - 0.6], color=DARK, lw=0.8, zorder=6)
text(gx + 1.8, gy, "not heritable: stop", color=DARK, size=FS_SMALL, ha="left",
     fontweight="bold")

# ------------------------------------------------------------------ 4. forward selection
steps = DATA["forward_selection"]["steps"]
sx0, sx1 = 67.6, 92.4
sy0, sy1 = YL, TIP_Y[0] - 1.0  # baseline, top of axis
RMAX = 0.75


def ry(r):
    return sy0 - (sy0 - sy1) * r / RMAX


line([sx0, sx0, sx1], [sy1, sy0, sy0], color=GREY, lw=0.5)
text(sx0 + 0.5, sy1, "R²", ha="left", va="top")
text(sx0, YL + 1.6, "\n".join(wrap(R2_DEF, sx1 - sx0, FS_DEF, "normal")), color=GREY,
     size=FS_DEF, ha="left", va="top", linespacing=1.15)
xs = [sx0 + 0.4, sx0 + 3.4, sx0 + 11.4, sx1 - 0.8]
r = [s["r2"] for s in steps]
line([xs[0], xs[1]], [ry(r[0])] * 2)
line([xs[1]] * 2, [ry(r[0]), ry(r[1])], color=CLOSE, lw=1.5)
line([xs[1], xs[2]], [ry(r[1])] * 2)
line([xs[2]] * 2, [ry(r[1]), ry(r[2])], color=OPEN, lw=1.5)
line([xs[2], xs[3]], [ry(r[2])] * 2, color=GREY, linestyle=(0, (1.5, 1)))
text(xs[3], ry(r[2]), "×", color=GREY, size=FS_SUB)


def pct(d):
    return f"+{100 * d:.1f}%" if d >= 0.01 else f"+{100 * d:.2f}%"


ym = ry(r[1] / 2)
text(xs[1] + 0.6, ym - 1.1, f"+{LABEL['RE112']}", color=CLOSE, ha="left")
text(xs[1] + 0.6, ym + 1.1, f"ΔR² {pct(steps[1]['delta_r2'])}", color=CLOSE,
     ha="left", size=FS_SMALL)
ym = ry(r[1]) + 2.2
text(xs[2] + 0.6, ym - 1.1, f"+{LABEL['RE63']}", color=OPEN, ha="left")
text(xs[2] + 0.6, ym + 1.1, f"ΔR² {pct(steps[2]['delta_r2'])}", color=OPEN,
     ha="left", size=FS_SMALL)
text(xs[3] + 0.3, ry(r[2]) - 1.5, f"next {pct(steps[3]['delta_r2'])}: stop",
     color=GREY, size=FS_SMALL, ha="right")

# ------------------------------------------------------------------ 5. drivers
tx = draw_tree(96.4, 6.8, drivers=True)
for i, t in enumerate(TIPS):
    cell(tx + 0.6, TIP_Y[i], t["actuation"])
LX = tx + 0.6 + CELL_W + 0.7
drv = {d["site"]: d for d in DATA["joint_model"]["drivers"]}


def p_fmt(p):
    m, e = f"{p:.0e}".split("e")
    return f"P = {m}e−{int(e[1:])}"


for site, c, word, shift in (("RE63", OPEN, "opening", -0.9), ("RE112", CLOSE, "closing", 0.4)):
    d = drv[site]
    ys = [TIP_Y[i] for i, t in enumerate(TIPS) if t[site]]
    yc = (min(ys) + max(ys)) / 2 + shift
    sign = "+" if d["beta"] > 0 else "−"
    rows = [(LABEL[site], FS_SUB, "bold"), (f"{d['change']}, {word}", FS_SMALL, "normal"),
            (f"β = {sign}{abs(d['beta']):.2f}", FS_SMALL, "normal"),
            (p_fmt(d["p"]), FS_SMALL, "normal")]
    for j, (s, fs, fw) in enumerate(rows):
        text(LX, yc + (j - 1.5) * 1.95, s, size=fs, color=c, ha="left", fontweight=fw)

text(W - 0.6, H - 1.0, f"example: {DATA['peak']['id'].split('_')[0]}:"
     f"{int(DATA['peak']['id'].split('_')[1]):,}", size=5.0, color=GREY, ha="right")

for x in (21.0, 45.4, gx, 94.4):
    chevron(x)

for ext, kw in (("svg", {}), ("pdf", {}), ("png", {"dpi": 300})):
    fig.savefig(f"{OUT}.{ext}", facecolor="white", **kw)
print(f"wrote {OUT}.svg/.pdf/.png ({W:.0f} x {H:.1f} mm); title lines:",
      [len(t) for t in TITLE_LINES])
