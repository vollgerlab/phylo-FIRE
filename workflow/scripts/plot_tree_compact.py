#!/usr/bin/env python3
"""Compact figure panel: tree | % FIRE | driver carriers | per-driver densities.

A main-figure-sized version of plot_tree_actuation_html.py without the MSA.
Driver selection, ordering (by p_fdr), colours and RE-relative labels match
the full plot, so the two can sit side by side.

Usage (paths from a phylo-FIRE run, e.g. results/test):
    pixi run python workflow/scripts/plot_tree_compact.py \\
        --tree      results/test/trees/chr8_70445067_70445271.rooted.nwk \\
        --actuation results/test/fastas/chr8_70445067_70445271.tsv \\
        --phylo-mm  results/test/phylo_mm/chr8_70445067_70445271.tsv \\
        --alignment results/test/aligned/chr8_70445067_70445271.fa --flank 250 \\
        --out       chr8_70445067_70445271.compact.pdf
"""

import argparse
import csv
import gzip
from pathlib import Path

import matplotlib
import matplotlib.patches

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from Bio import Phylo  # noqa: E402

from plot_tree_actuation_html import (  # noqa: E402
    ACT_COLORS,
    ACT_MISSING,
    SIG_GROUP_PALETTE,
    _branch_is_sig,
    compute_positions,
    load_actuation,
    load_alignment,
    load_phylo_mm,
    parse_leaf_meta,
    re_band_columns,
)

MARKERS = ["o"]  # colour distinguishes drivers; keep every marker a circle


def act_color(a):
    """10-bin actuation colour (same bins as the full plot)."""
    return ACT_MISSING if a is None else ACT_COLORS[min(int(a * 10), 9)]


def kde(x, grid):
    """Gaussian KDE with Silverman's bandwidth, peak-normalised to 1."""
    x = np.asarray(x, float)
    if len(x) < 2 or np.std(x) == 0:
        y = np.exp(-0.5 * ((grid - x.mean()) / 2.0) ** 2)
    else:
        bw = 1.06 * np.std(x, ddof=1) * len(x) ** -0.2
        y = np.exp(-0.5 * ((grid[:, None] - x[None, :]) / bw) ** 2).sum(1)
    return y / y.max()


def collapse_zero(clade, eps=1e-5):
    """Splice out internal children with branch length < eps (bottom-up), making polytomies."""
    for c in clade.clades:
        collapse_zero(c, eps)
    kids = []
    for c in clade.clades:
        if not c.is_terminal() and (c.branch_length or 0) < eps:
            for g in c.clades:
                g.branch_length = (g.branch_length or 0) + (c.branch_length or 0)
                kids.append(g)
        else:
            kids.append(c)
    clade.clades = kids


def drivers(pmm, tips, alignment, flank):
    """Significant driver branches -> [{label, color, carriers, beta, p}], ordered by p_fdr."""
    sig = sorted(
        (b for b in pmm["branches"] if _branch_is_sig(b)),
        key=lambda b: b["p_fdr"] if isinstance(b.get("p_fdr"), float) else 1.0,
    )
    band = None
    if alignment:
        aln_len = max(len(s) for s in alignment.values())
        band = re_band_columns(alignment, tips, flank, aln_len)
    out = []
    for i, b in enumerate(sig):
        sites = [int(s) for s in b["linked_sites"] if s.isdigit()]
        if sites and band:
            label = f"RE{min(sites) - band[0]}"
        elif sites:
            label = f"col{min(sites)}"
        else:
            label = f"d{i}"
        out.append({
            "label": label,
            "color": SIG_GROUP_PALETTE[i % len(SIG_GROUP_PALETTE)],
            "marker": MARKERS[i % len(MARKERS)],
            "carriers": {parse_leaf_meta(d)["sample_hap"] for d in b["descendants"]},
            "beta": b.get("beta"),
            "p": b.get("p") if isinstance(b.get("p"), float) else b.get("p_fdr"),
            "p_fdr": b.get("p_fdr"),
        })
    # Nested driver -> compare against the smallest enclosing driver's other carriers
    # (the full plot's local contrast); otherwise against everyone else.
    for d in out:
        parents = [o for o in out if o is not d and o["carriers"] > d["carriers"]]
        d["parent"] = min(parents, key=lambda o: len(o["carriers"]))["carriers"] if parents else None
    return out


def load_variants(vcf: Path, peak: str) -> dict:
    """RE position -> {hg38, ref, alt, rsid} for one peak from variants.hg38[.rsid].vcf.gz.

    ALT is the driver allele when the record has DALLELE, else the VCF ALT.
    """
    out = {}
    with gzip.open(vcf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            info = dict(kv.split("=", 1) if "=" in kv else (kv, True) for kv in f[7].split(";"))
            if info.get("PEAK") != peak or "RE" not in info:
                continue
            out[f"RE{info['RE']}"] = {
                "hg38": f"{f[0]}:{int(f[1]):,}", "ref": f[3], "alt": info.get("DALLELE", f[4]),
                "rsid": f[2] if f[2] not in (".", "") else "",
            }
    return out


def write_source_data(out: Path, actuation: Path, sh: list, drv: list, meta: dict) -> None:
    """Write the numbers behind the figure next to it.

    <stem>.source_data.tsv: one row per haplotype in plot (top-to-bottom) order.
      For each driver: carrier (1/0) and group = the density side it is drawn on
      ("carrier", "comparison", or "" when excluded from that driver's contrast).
    <stem>.drivers.tsv: one row per driver with the values printed on the figure.
    """
    with actuation.open() as fh:
        rows = {r["sample"]: r for r in csv.DictReader(fh, delimiter="\t")}
    stem = out.with_suffix("")
    peak = meta.get("peak_id", "")

    def num(x):
        return "" if x in (None, "", "NA") else x

    with open(f"{stem}.source_data.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        head = ["peak_id", "plot_order", "haplotype", "fire_coverage", "coverage", "pct_fire"]
        for d in drv:
            head += [f"{d['label']}_carrier", f"{d['label']}_group"]
        w.writerow(head)
        for i, s in enumerate(sh, 1):
            r = rows.get(s, {})
            a = r.get("actuation", "NA")
            pct = "" if a in ("", "NA") else f"{float(a) * 100:.2f}"
            line = [peak, i, s, num(r.get("fire_coverage")), num(r.get("coverage")), pct]
            for d in drv:
                grp = "carrier" if s in d["carriers"] else ("comparison" if s in d.get("comp", ()) else "")
                line += [int(s in d["carriers"]), grp]
            w.writerow(line)

    act = {s: float(r["actuation"]) * 100 for s, r in rows.items() if r.get("actuation") not in ("", "NA", None)}
    with open(f"{stem}.drivers.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["peak_id", "driver", "n_carriers", "n_comparison", "mean_pct_fire_carriers",
                    "mean_pct_fire_comparison", "delta_pct_fire", "beta", "p", "p_fdr",
                    "hg38_pos", "ref", "driver_allele", "rsid"])
        for d in drv:
            car = [act[s] for s in d["carriers"] if s in act]
            comp = [act[s] for s in d.get("comp", ()) if s in act]
            mc = np.mean(car) if car else float("nan")
            mr = np.mean(comp) if comp else float("nan")
            w.writerow([peak, d["label"], len(car), len(comp), f"{mc:.2f}", f"{mr:.2f}", f"{mc - mr:.2f}",
                        num(d["beta"]), num(d["p"]), num(d.get("p_fdr")),
                        *(d.get("var", {}).get(k, "") for k in ("hg38", "ref", "alt", "rsid"))])
    print(f"wrote {stem}.source_data.tsv, {stem}.drivers.tsv")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tree", required=True, type=Path)
    ap.add_argument("--actuation", required=True, type=Path)
    ap.add_argument("--phylo-mm", required=True, type=Path)
    ap.add_argument("--alignment", type=Path, help="flanked MSA; only used to label drivers RE<pos>")
    ap.add_argument("--flank", type=int, default=0)
    ap.add_argument("--variants", type=Path,
                    help="variants.hg38.rsid.vcf.gz (or variants.hg38.vcf.gz) to label drivers with change + rsID")
    ap.add_argument("--out", required=True, type=Path, help=".pdf/.png/.svg")
    ap.add_argument("--title", default=None, help="default: peak id from phylo_mm header")
    ap.add_argument("--no-labels", action="store_true", help="drop the (tiny, editable) haplotype names")
    ap.add_argument("--width", type=float, default=3.0, help="figure width, inches")
    ap.add_argument("--height", type=float, default=3.0, help="figure height, inches")
    args = ap.parse_args()

    tree = Phylo.read(str(args.tree), "newick")
    # IQ-TREE resolves identical sequences into an arbitrary chain of zero-length
    # branches, which scatters near-identical singletons (e.g. one flank SNV) away
    # from their group. Collapse those into polytomies, then order each node's
    # children: lone leaves first (shortest branch first), then clades by size.
    collapse_zero(tree.root)
    act, _ = load_actuation(args.actuation)

    def leaf_act(c):
        a = act.get(parse_leaf_meta(c.name)["sample_hap"], act.get(c.name))
        return -1.0 if a is None else a

    # siblings: identical leaves, then leaves on their own branch, then sub-clades
    # by size; leaves within each group sorted by % FIRE (highest on top)
    for clade in tree.find_clades(order="postorder"):
        clade.clades.sort(key=lambda c: (
            (0 if (c.branch_length or 0) < 1e-5 else 1, -leaf_act(c)) if c.is_terminal()
            else (2, c.count_terminals())))
    pos = compute_positions(tree.root)
    tips = [t.name for t in tree.get_terminals()]
    n = len(tips)
    sh = [parse_leaf_meta(t)["sample_hap"] for t in tips]
    acts = [act.get(s, act.get(t)) for s, t in zip(sh, tips)]
    pmm = load_phylo_mm(args.phylo_mm)
    aln = load_alignment(args.alignment) if args.alignment else None
    drv = drivers(pmm, tips, aln, args.flank)
    if args.variants:
        peak_id = str(pmm["meta"].get("peak_id", ""))
        var = load_variants(args.variants, peak_id)
        for d in drv:
            d["var"] = var.get(d["label"], {})
    nd = max(len(drv), 1)

    # ── layout in inches: [tree][labels?][%FIRE][markers]  gap  [densities / colorbar]
    W, H = args.width, args.height
    fig = plt.figure(figsize=(W, H))
    # Arial is on macOS and in Illustrator, so embedded/linked PDFs keep their fonts;
    # Type 42 keeps text editable.
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                         "pdf.fonttype": 42, "ps.fonttype": 42, "font.size": 5, "axes.linewidth": 0.4,
                         "xtick.major.width": 0.4, "xtick.major.size": 1.5, "xtick.major.pad": 1})

    def ax_in(x, y, w, h):
        return fig.add_axes([x / W, y / H, w / W, h / H])

    bot, top = 0.34, H - 0.32               # leaf-row band
    tree_w, lab_w, act_w, mk_w = 0.42 * W, (0.0 if args.no_labels else 0.2), 0.10, 0.07 * nd
    x = 0.05
    ax_tree = ax_in(x, bot, tree_w, top - bot); x += tree_w
    ax_lab = ax_in(x, bot, lab_w or 1e-3, top - bot); x += lab_w
    ax_act = ax_in(x + 0.02, bot, act_w, top - bot); x += act_w + 0.02
    ax_mk = ax_in(x + 0.01, bot, mk_w, top - bot); x += mk_w + 0.01
    den_x = x + 0.10
    den_w = W - den_x - 0.10  # densities get whatever the tree leaves
    den_h = min(top - bot - 0.45, 0.75 * nd)
    ax_den = ax_in(den_x, top - den_h, den_w, den_h)
    for ax in (ax_tree, ax_lab, ax_act, ax_mk):
        ax.set_ylim(n - 0.5, -0.5)
        ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)

    # tree: colour each driver's stem + its carriers' tip extensions
    xmax = max(p["x"] for p in pos.values()) or 1.0
    tip_color, stem_color = {}, {}
    for d in reversed(drv):  # most significant driver wins overlaps
        car = [t for t in tree.get_terminals() if parse_leaf_meta(t.name)["sample_hap"] in d["carriers"]]
        for t in car:
            tip_color[t.name] = d["color"]
        stem_color[id(tree.common_ancestor(car))] = d["color"]
    for clade in tree.find_clades():
        p = pos[id(clade)]
        for c in clade.clades:
            q = pos[id(c)]
            ax_tree.plot([p["x"], q["x"]], [q["y"]] * 2, color=stem_color.get(id(c), "#444444"),
                         lw=1.6 if id(c) in stem_color else 0.4, solid_capstyle="butt")
        if clade.clades:
            ys = [pos[id(c)]["y"] for c in clade.clades]
            ax_tree.plot([p["x"]] * 2, [min(ys), max(ys)], color="#444444", lw=0.4)
        if clade.is_terminal():
            ax_tree.plot([p["x"], xmax * 1.02], [p["y"]] * 2, ls=":", lw=0.3,
                         color=tip_color.get(clade.name, "#d0d0d0"))
    ax_tree.set_xlim(0, xmax * 1.02)
    ax_tree.set_xticks([0, xmax])
    ax_tree.set_xticklabels(["0", f"{xmax:.3f}"])
    ax_tree.set_xlabel("branch length", labelpad=1)
    ax_tree.spines["bottom"].set_visible(True)

    if not args.no_labels:
        fs = min(5.0, 72 * (top - bot) / n * 0.8)  # ~1 row pitch; tiny at 3 in
        for i, s in enumerate(sh):
            ax_lab.text(0.02, i, s, va="center", ha="left", fontsize=fs, color=tip_color.get(tips[i], "#222222"))
    ax_lab.set_xlim(0, 1)
    ax_lab.set_xticks([])

    ax_act.barh(range(n), [1] * n, height=1.0, color=[act_color(a) for a in acts], lw=0)
    ax_act.set_xlim(0, 1)
    ax_act.set_xticks([])
    ax_act.set_xlabel("%\nFIRE", labelpad=1, linespacing=0.9)

    for j, d in enumerate(drv):
        ys = [i for i, s in enumerate(sh) if s in d["carriers"]]
        ax_mk.scatter([j] * len(ys), ys, marker=d["marker"], s=2.5, color=d["color"], lw=0)
    ax_mk.set_xlim(-0.6, nd - 0.4)
    ax_mk.set_xticks(range(len(drv)))
    ax_mk.set_xticklabels([d["label"] for d in drv], rotation=90, fontsize=4.5)
    ax_mk.tick_params(length=0)
    for lab, d in zip(ax_mk.get_xticklabels(), drv):
        lab.set_color(d["color"])

    # densities: one row per driver, carriers below the line, comparison above
    grid = np.linspace(0, 100, 400)
    sh_act = {s: a * 100 for s, a in zip(sh, acts) if a is not None}
    rng = np.random.default_rng(0)
    for j, d in enumerate(drv):
        comp = (d["parent"] - d["carriers"]) if d["parent"] else set(sh_act) - d["carriers"]
        d["comp"] = comp
        car = [sh_act[s] for s in d["carriers"] if s in sh_act]
        rest = [sh_act[s] for s in comp if s in sh_act]
        y0 = -j * 1.75
        if rest:
            ax_den.fill_between(grid, y0, y0 + 0.5 * kde(rest, grid), color="#bbbbbb", alpha=0.4, lw=0)
            ax_den.scatter(rest, y0 + 0.06 + 0.32 * rng.random(len(rest)), s=0.8,
                           facecolors="none", edgecolors="#888888", lw=0.25)
            ax_den.plot([np.mean(rest)] * 2, [y0, y0 + 0.5], color="#666666", ls="--", lw=0.5)
        if car:
            ax_den.fill_between(grid, y0 - 0.5 * kde(car, grid), y0, color=d["color"], alpha=0.3, lw=0)
            ax_den.scatter(car, y0 - 0.06 - 0.32 * rng.random(len(car)), s=2, marker=d["marker"],
                           color=d["color"], lw=0)
            ax_den.plot([np.mean(car)] * 2, [y0 - 0.5, y0], color=d["color"], ls="--", lw=0.5)
        ax_den.axhline(y0, color="#999999", lw=0.3)
        v = d.get("var", {})
        change = f" {v['ref']}>{v['alt']}" if v.get("ref") else ""
        head_txt = d["label"] + change + (f"  Δ={np.mean(car) - np.mean(rest):+.0f}%" if car and rest else "")
        stat_txt = "  ".join(
            ([v["rsid"]] if v.get("rsid") else [])
            + ([f"β={d['beta']:+.2f}"] if isinstance(d["beta"], float) else [])
            + ([f"p={d['p']:.0e}"] if isinstance(d["p"], float) else []))
        ax_den.text(0, y0 + 0.86, head_txt, ha="left", va="bottom", color=d["color"], fontweight="bold")
        ax_den.text(0, y0 + 0.56, stat_txt, ha="left", va="bottom", color=d["color"], fontsize=4.5)
    ax_den.set_xlim(-4, 104)  # pad so points at 0%/100% are not clipped
    ax_den.spines["bottom"].set_bounds(0, 100)
    ax_den.set_ylim(-(nd - 1) * 1.75 - 0.55, 1.2)
    ax_den.set_yticks([])
    ax_den.set_xticks([0, 50, 100])
    for s in ("left", "right", "top"):
        ax_den.spines[s].set_visible(False)
    ax_den.set_xlabel("% FIRE", labelpad=1)

    # %FIRE colour key for the column, under the densities
    ax_cb = ax_in(den_x, bot + 0.02, den_w, 0.07)
    for i, c in enumerate(ACT_COLORS):  # vector bins (an imshow raster breaks on Illustrator embed)
        ax_cb.add_patch(matplotlib.patches.Rectangle((i * 10, 0), 10, 1, fc=c, ec="none"))
    ax_cb.set_ylim(0, 1)
    ax_cb.set_xlim(-4, 104)  # same frame as the densities above
    ax_cb.add_patch(matplotlib.patches.Rectangle((0, 0), 100, 1, fill=False, lw=0.4, ec="#444444"))
    for s_ in ("left", "right", "top", "bottom"):
        ax_cb.spines[s_].set_visible(False)
    ax_cb.set_yticks([])
    ax_cb.set_xticks([0, 50, 100])
    ax_cb.set_xlabel("% FIRE (column)", labelpad=1)
    ax_cb.tick_params(length=1.5)
    # key drawn with real markers, not ○/● glyphs (missing from most fonts)
    ax_cb.scatter([1.5], [1.75], s=2, facecolors="none", edgecolors="#888888", lw=0.35, clip_on=False)
    ax_cb.text(4, 1.75, "grey = rest", va="center", fontsize=4, color="#555555", clip_on=False)
    ax_cb.scatter([40], [1.75], s=2.5, color="#555555", lw=0, clip_on=False)
    ax_cb.text(42.5, 1.75, "colour = carriers", va="center", fontsize=4, color="#555555", clip_on=False)

    meta = pmm["meta"]
    title = args.title or str(meta.get("peak_id", args.tree.stem.split(".")[0]))
    sub = [f"n={n}"]
    if isinstance(meta.get("n_unique_seqs"), float):
        sub.append(f"unique={int(meta['n_unique_seqs'])}")
    if isinstance(meta.get("h2_re_mean"), float):
        sub.append(f"H²(RE)={meta['h2_re_mean']:.2f}")
    fig.text(0.05 / W, 1 - 0.06 / H, title, ha="left", va="top", fontsize=6.5, fontweight="bold")
    fig.text(0.05 / W, 1 - 0.19 / H, " · ".join(sub), ha="left", va="top", fontsize=5, color="#444444")

    fig.savefig(args.out, dpi=600)
    print(f"wrote {args.out}")
    write_source_data(args.out, args.actuation, sh, drv, pmm["meta"])


if __name__ == "__main__":
    main()
