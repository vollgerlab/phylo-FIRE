#!/usr/bin/env python3
"""Plot a phylogenetic tree alongside a FIRE % actuation heatmap.

Layout (left to right):
  - Cladogram (IQ-TREE Newick, midpoint-rooted, ladderized)
  - One colored box per leaf: white (0%) -> dark red (100%), binned in 10% intervals
  - Legend

Usage:
    python plot_tree_actuation.py \\
        --tree trees/chr8_8453862_8454062.treefile \\
        --actuation fastas/chr8_8453862_8454062.tsv \\
        --output plots/chr8_8453862_8454062.pdf
"""

import argparse
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from Bio import Phylo

# ColorBrewer "Reds" 10-bin palette: white (0%) to dark red (100%)
BIN_COLORS = [
    "#fff5f0",  # 0-10%
    "#fee0d2",  # 10-20%
    "#fcbba1",  # 20-30%
    "#fc9272",  # 30-40%
    "#fb6a4a",  # 40-50%
    "#ef3b2c",  # 50-60%
    "#cb181d",  # 60-70%
    "#a50f15",  # 70-80%
    "#67000d",  # 80-90%
    "#3d0005",  # 90-100%
]
MISSING_COLOR = "#d0d0d0"
BIN_LABELS = [f"{i*10}-{i*10+10}%" for i in range(10)]


def actuation_color(value: float | None) -> str:
    if value is None:
        return MISSING_COLOR
    bin_idx = min(int(float(value) * 10), 9)
    return BIN_COLORS[bin_idx]


def load_actuation(tsv_path: Path) -> dict[str, float | None]:
    data: dict[str, float | None] = {}
    with tsv_path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            val = row.get("actuation", "NA")
            data[row["sample"]] = None if val in ("NA", "") else float(val)
    return data


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot phylogenetic tree with FIRE actuation heatmap.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--tree", required=True, type=Path, help="Newick tree file")
    parser.add_argument("--actuation", required=True, type=Path, help="Per-sample actuation TSV")
    parser.add_argument("--output", required=True, type=Path, help="Output PDF or PNG path")
    parser.add_argument("--peak-id", type=str, default=None, help="Peak ID for plot title")
    args = parser.parse_args()

    tree = Phylo.read(str(args.tree), "newick")
    tree.root_at_midpoint()
    tree.ladderize()

    actuation = load_actuation(args.actuation)

    # Bio.Phylo assigns y = 1..n to terminals in the order of tree.get_terminals()
    leaves = [c.name for c in tree.get_terminals()]
    n = len(leaves)

    # Size figure to the number of leaves; keep a minimum height for legends
    fig_height = max(5, n * 0.22)
    fig, (ax_tree, ax_heat, ax_leg) = plt.subplots(
        1, 3,
        figsize=(14, fig_height),
        gridspec_kw={"width_ratios": [5, 0.5, 1.5], "wspace": 0.02},
    )

    # Draw tree
    Phylo.draw(
        tree,
        axes=ax_tree,
        do_show=False,
        label_func=lambda c: c.name if c.is_terminal() else "",
    )
    ax_tree.set_ylabel("")
    ax_tree.set_xlabel("branch length", fontsize=8)
    title = args.peak_id or args.tree.stem
    ax_tree.set_title(title, fontsize=9, pad=6)
    ax_tree.tick_params(axis="y", labelsize=6)

    # Sync y-limits so heatmap rows align with tree leaves
    ylim = ax_tree.get_ylim()
    ax_heat.set_ylim(ylim)

    # Draw one colored box per leaf
    # Bio.Phylo places leaf i (0-indexed from get_terminals()) at y = i+1,
    # then inverts so y=1 appears at the top of the figure.
    ax_heat.set_xlim(0, 1)
    ax_heat.axis("off")
    ax_heat.set_title("% act.", fontsize=7, pad=4)

    box_h = 0.75
    for i, sample in enumerate(leaves):
        y = i + 1
        val = actuation.get(sample, None)
        color = actuation_color(val)
        rect = mpatches.FancyBboxPatch(
            (0.1, y - box_h / 2),
            0.8,
            box_h,
            boxstyle="square,pad=0",
            facecolor=color,
            edgecolor="#aaaaaa",
            linewidth=0.4,
        )
        ax_heat.add_patch(rect)

    # Legend
    ax_leg.axis("off")
    legend_patches = [
        mpatches.Patch(
            facecolor=BIN_COLORS[i],
            edgecolor="#aaaaaa",
            linewidth=0.5,
            label=BIN_LABELS[i],
        )
        for i in range(10)
    ]
    legend_patches.append(
        mpatches.Patch(
            facecolor=MISSING_COLOR,
            edgecolor="#aaaaaa",
            linewidth=0.5,
            label="missing",
        )
    )
    ax_leg.legend(
        handles=legend_patches,
        loc="center left",
        frameon=False,
        fontsize=7,
        title="FIRE actuation",
        title_fontsize=7.5,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(args.output), bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"Saved {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
