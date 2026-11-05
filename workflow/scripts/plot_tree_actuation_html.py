#!/usr/bin/env python3
"""Interactive phylogenetic tree + FIRE actuation + MSA visualization via Plotly.

Layout (left to right):
  tree | leaf labels | PCLAI ancestry stripe | FIRE % actuation | MSA
                                              with sequence logo above

Writes a PDF via kaleido, plus interactive HTML with --html.

The PCLAI stripe (lookup in pclai.py) uses HPRC R2's pre-baked palette. It is
omitted when --pclai-bed is not supplied or no leaf matches (e.g. peaks whose
chrom is not a CHM13 chromosome).

Usage:
    python plot_tree_actuation_html.py \\
        --tree      trees/chr8_8453862_8454062.treefile \\
        --actuation fastas/chr8_8453862_8454062.tsv \\
        --alignment aligned/chr8_8453862_8454062.fa \\
        --html      plots/chr8_8453862_8454062.html \\
        --pdf       plots/chr8_8453862_8454062.pdf \\
        --pclai-bed resources/pclai/hprc_r2.bed.gz   # optional
"""

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path

from Bio import Phylo
import plotly.graph_objects as go

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common_io import read_fasta  # noqa: E402
from pclai import PCLAI_MISSING, PCLAI_NEAREST_BP, query_pclai  # noqa: E402

# ── Colour palettes ───────────────────────────────────────────────────────────

# Actuation 10-bin palette: pure white at 0-10%, then ColorBrewer Reds[9].
# Pure-white bottom bin makes 0-10% bars effectively blank against the white
# background — the outline still draws so the row stays locatable but the
# colour signal only appears once actuation crosses 10%.
ACT_COLORS = [
    "#ffffff", "#fff5f0", "#fee0d2", "#fcbba1", "#fc9272",
    "#fb6a4a", "#ef3b2c", "#cb181d", "#a50f15", "#67000d",
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
MSA_INVARIANT = "#EEEEEE"   # near-white — base matches column consensus
MSA_MISSING   = "#555555"   # very dark grey — gap / N / unmapped

# MSA z-encoding:  0=var-A  1=var-T  2=var-G  3=var-C  4=invariant  5=missing
# Colorscale: one entry per integer z value, normalized to [0,1] via zmax=5.
# Must start at exactly 0.0 and end at exactly 1.0 — Plotly silently rejects
# colorscales that don't span [0,1] and falls back to "Reds".
BASE_TO_CODE = {"A": 0, "T": 1, "G": 2, "C": 3}
DNA_COLORSCALE = [
    [0.0, "#009600"],       # z=0: variant A — green
    [0.2, "#FF4444"],       # z=1: variant T — red
    [0.4, "#FFB300"],       # z=2: variant G — amber
    [0.6, "#4444FF"],       # z=3: variant C — blue
    [0.8, MSA_INVARIANT],   # z=4: matches consensus — near-white
    [1.0, MSA_MISSING],     # z=5: gap / N — very dark grey
]

# Categorical palette for paralog groups (ColorBrewer Set1 + Set2 mix).
PARALOG_PALETTE = [
    "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
    "#a65628", "#f781bf", "#66c2a5", "#fc8d62", "#8da0cb",
]
LEAF_DEFAULT = "#222222"
INTERNAL_BRANCH = "#666666"


# ── Leaf metadata parsing ─────────────────────────────────────────────────────

def parse_leaf_meta(name: str) -> dict:
    """Parse leaf names into (locus, paralog, sample_hap).

    Three formats are supported:
      - paralog peaks:    "LOCUS__PARALOG__SAMPLE_HAP" (PARE)
      - tagged union:     "TAG__SAMPLE_HAP"  (e.g. extract_fastas --highlight-file)
      - plain union:      "SAMPLE_HAP"
    The ``paralog`` field is what drives leaf-label colour + the paralog
    stripe, so the 2-segment form is the natural way to highlight a subset
    of leaves without dragging in a per-leaf locus token.
    """
    parts = name.split("__")
    if len(parts) == 3:
        return {"locus": parts[0], "paralog": parts[1], "sample_hap": parts[2]}
    if len(parts) == 2:
        return {"locus": None, "paralog": parts[0], "sample_hap": parts[1]}
    return {"locus": None, "paralog": None, "sample_hap": name}


def assign_paralog_colors(metas: list[dict]) -> tuple[dict[str, str], list[str], bool]:
    """Map paralog → stable colour and return per-leaf label colours.

    Returns (paralog_color_map, per_leaf_colors, has_paralogs).
    """
    paralogs = sorted({m["paralog"] for m in metas if m["paralog"]})
    color_map = {p: PARALOG_PALETTE[i % len(PARALOG_PALETTE)] for i, p in enumerate(paralogs)}
    leaf_colors = [color_map.get(m["paralog"], LEAF_DEFAULT) for m in metas]
    return color_map, leaf_colors, bool(paralogs)


# ── Data loaders ──────────────────────────────────────────────────────────────

def load_actuation(path: Path) -> tuple[dict[str, float | None], dict[str, int | None]]:
    """Return (actuation, coverage) dicts keyed by sample name."""
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


def load_phylo_mm(path: Path) -> dict:
    """Parse a phylo_mm TSV into ``{"meta": header k=v, "branches": per-variant rows}``.

    Rows carry the marginal carriers-vs-rest stats (mean_in/out, t_stat, p_raw,
    p_fdr), the joint forward-selection stats (beta, se, ci_*, p) and
    ``is_split`` (selected driver). Missing or malformed file → empty.
    """
    out: dict = {"meta": {}, "branches": []}
    try:
        with path.open() as fh:
            lines = [ln.rstrip("\n") for ln in fh]
    except (OSError, AttributeError):
        return out
    data_lines: list[str] = []
    for ln in lines:
        if not ln.strip():
            continue
        if ln.startswith("#"):
            kv = ln.lstrip("# ").split("=", 1)
            if len(kv) == 2:
                key, val = kv[0].strip(), kv[1].strip()
                try:
                    out["meta"][key] = float(val)
                except ValueError:
                    out["meta"][key] = val
            continue
        data_lines.append(ln)
    if len(data_lines) < 2:
        return out
    header, *rows = data_lines
    cols = header.split("\t")
    for raw in rows:
        f = raw.split("\t")
        if len(f) != len(cols):
            continue
        rec = dict(zip(cols, f))
        for k in ("branch_id", "n_below"):
            try:
                rec[k] = int(rec[k])
            except (KeyError, ValueError):
                pass
        for k in ("mean_in", "mean_out", "t_stat", "p_raw", "p_fdr",
                  "lasso_coef", "beta", "se",
                  "ci_lower", "ci_upper", "p"):
            raw_val = rec.get(k, "")
            if raw_val in ("", "NA"):
                rec[k] = None
            else:
                try:
                    rec[k] = float(raw_val)
                except ValueError:
                    rec[k] = None
        rec["selected_by_lasso"] = (
            str(rec.get("selected_by_lasso", "")).upper() == "TRUE"
        )
        rec["is_split"] = str(rec.get("is_split", "")).upper() == "TRUE"
        rec["descendants"] = tuple(
            d for d in rec.get("descendants", "").split(",") if d
        )
        rec["linked_sites"] = tuple(
            s for s in rec.get("linked_sites", "").split(",") if s
        )
        out["branches"].append(rec)
    return out


def load_alignment(path: Path) -> dict[str, str]:
    """Read a FASTA into {name: sequence}, stripping MAFFT's ``_R_`` prefix
    (``--adjustdirection``) so names match the actuation TSV."""
    return read_fasta(path, strip_r=True)


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


def branch_xy_split(root, pos) -> tuple[list, list, dict[int, tuple[list, list]]]:
    """Split branches into internal segments + per-leaf terminal segments.

    Returns:
      internal_xs, internal_ys: vertical clade-joiners + horizontal internal branches
      terminals: dict[id(leaf_clade)] -> (xs, ys) for that leaf's horizontal branch
    """
    internal_xs: list = []
    internal_ys: list = []
    terminals: dict[int, tuple[list, list]] = {}
    for clade in root.find_clades(order="level"):
        if not clade.clades:
            continue
        px = pos[id(clade)]["x"]
        child_ys = [pos[id(c)]["y"] for c in clade.clades]
        # Vertical bar joining children (always internal)
        internal_xs += [px, px, None]
        internal_ys += [min(child_ys), max(child_ys), None]
        # Horizontal branch to each child: terminal if child is a leaf
        for child in clade.clades:
            seg_xs = [px, pos[id(child)]["x"]]
            seg_ys = [pos[id(child)]["y"], pos[id(child)]["y"]]
            if child.is_terminal():
                terminals[id(child)] = (seg_xs, seg_ys)
            else:
                internal_xs += seg_xs + [None]
                internal_ys += seg_ys + [None]
    return internal_xs, internal_ys, terminals


# ── Sequence logo helpers ─────────────────────────────────────────────────────

def re_band_columns(
    alignment: dict[str, str], names: list[str], flank: int, aln_len: int
) -> tuple[int, int] | None:
    """Alignment-column span [lo, hi] of the regulatory element.

    When `flank` bp of context was added each side (rule add_flanks), the RE
    occupies the interior of every sequence: skip `flank` non-gap characters
    from each end, and the remaining non-gap span is the RE. The union across
    rows (min left column, max right column) is the RE band. Returns None when
    flanking is off or no row has enough sequence to locate the RE.
    """
    if flank <= 0:
        return None
    lo, hi = aln_len, -1
    for nm in names:
        s = alignment.get(nm)
        if not s:
            continue
        nongap = [i for i, c in enumerate(s) if c != "-"]
        if len(nongap) <= 2 * flank:
            continue   # flanks consumed the whole sequence — can't locate RE
        left_col = nongap[flank]
        right_col = nongap[len(nongap) - flank - 1]
        lo = min(lo, left_col)
        hi = max(hi, right_col)
    return (lo, hi) if hi >= lo else None


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


# ── Drivers ───────────────────────────────────────────────────────────────────

# Driver palette (Plotly D3 qualitative). Red is the first driver (smallest
# marginal p_fdr); later drivers cycle through the rest.
SIG_GROUP_PALETTE = [
    "#d62728", "#1f77b4", "#2ca02c", "#9467bd", "#ff7f0e",
    "#17becf", "#e377c2", "#8c564b", "#bcbd22", "#7f7f7f",
]


def _branch_is_sig(b: dict) -> bool:
    """Draw a variant iff it's a contributing variant (`is_split`).

    `is_split` marks the variants the forward-selection joint regression kept —
    those that explain actuation *beyond* the variants already in the model
    (joint p < 0.05). We gate on that alone, NOT on the marginal p_fdr: a nested
    variant that switches actuation back off within a stronger variant's clade
    has a real joint contribution but a muddied marginal carriers-vs-rest effect.
    """
    return bool(b.get("is_split"))


DRIVER_SYMBOLS = ["cross", "circle", "star", "diamond", "triangle-up",
                  "x", "square", "hexagram"]


def driver_groups(phylo_mm: dict | None) -> list[dict]:
    """Drivers (is_split rows) ordered by marginal p_fdr, smallest first.

    Each entry is ``{"rec", "idx", "color", "symbol", "links"}``: the phylo_mm
    row, its position among the is_split rows, its palette colour (index 0 =
    red), its marker symbol and its linked MSA sites. phylo_mm merges identical
    and complementary partitions, so every driver has its own carrier set.
    """
    sig = [b for b in (phylo_mm or {}).get("branches") or [] if _branch_is_sig(b)]
    ranked = sorted(
        enumerate(sig),
        key=lambda ib: ib[1]["p_fdr"] if isinstance(ib[1].get("p_fdr"), float) else 1.0,
    )
    return [
        {
            "rec": b,
            "idx": idx,
            "color": SIG_GROUP_PALETTE[i % len(SIG_GROUP_PALETTE)],
            "symbol": DRIVER_SYMBOLS[i % len(DRIVER_SYMBOLS)],
            "links": frozenset(b.get("linked_sites") or ()),
        }
        for i, (idx, b) in enumerate(ranked)
    ]


# Synthetic mutagenesis constructs (from make_mutagenesis_panel --aligned-out)
# render as grey pseudo-leaf rows for visual validation: same MSA column frame,
# no tree branch, no actuation. A distinct label colour sets them apart.
PANEL_LABEL_COLOR = "#1f6f8b"


class _PanelRow:
    """A synthetic construct shown as a pseudo-leaf row next to its parent."""
    __slots__ = ("name", "label")

    def __init__(self, name: str, label: str):
        self.name = name
        self.label = label


def load_panel(path: Path) -> list[dict]:
    """make_mutagenesis_panel --aligned-out FASTA → [{name, seq, anchor, label}].

    Each header is ``>NAME anchor_leaf=… type=… edit=… label=…``; the sequence
    is in the SAME alignment column frame as the main MSA (so it slots in)."""
    rows: list[dict] = []
    name, hdr, buf = None, {}, []

    def add():
        rows.append({
            "name": name, "seq": "".join(buf),
            "anchor": hdr.get("anchor_leaf", ""), "label": hdr.get("label", name),
        })

    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                add()
            parts = line[1:].split()
            name = parts[0]
            hdr = dict(p.split("=", 1) for p in parts[1:] if "=" in p)
            buf = []
        else:
            buf.append(line)
    if name is not None:
        add()
    return rows

# ── Formatting helpers ────────────────────────────────────────────────────────

def _fmt_p(p) -> str:
    """Compact p-value (None → "NA")."""
    return f"{p:.2g}" if isinstance(p, float) else "NA"


def _fmt_beta(b) -> str:
    """Signed effect size (None → "NA")."""
    return f"{b:+.2f}" if isinstance(b, float) else "NA"


def _marginal_delta(rec: dict) -> float | None:
    """Marginal Δ = mean_in − mean_out: carriers vs all other haplotypes."""
    mi, mo = rec.get("mean_in"), rec.get("mean_out")
    if isinstance(mi, float) and isinstance(mo, float):
        return mi - mo
    return None


def _hex_rgba(h: str, a: float) -> str:
    """``#rrggbb`` → ``rgba(r,g,b,a)``; other strings pass through without '#'."""
    h = h.lstrip("#")
    if len(h) == 6:
        r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
        return f"rgba({r},{g},{b},{a})"
    return h


def _site_label(s, re_band: tuple[int, int] | None) -> str:
    """Label a linked site by RE-relative position when flanks are present.

    phylo_mm reports linked sites as 1-based alignment columns; subtracting the
    RE's 0-based start column gives the 1-based position within the RE (`RE13`
    = 13 bp into the element). Without flanks the column is shown as-is.
    """
    if re_band is None:
        return str(s)
    try:
        return f"RE{int(s) - re_band[0]}"
    except (TypeError, ValueError):
        return str(s)


def _sample_hap(name: str) -> str:
    return parse_leaf_meta(name)["sample_hap"]


# ── Row order ─────────────────────────────────────────────────────────────────

def prepare_tree(tree) -> None:
    """Ladderize and strip MAFFT's ``_R_`` prefix from leaf names, in place.

    The tree is already midpoint-rooted upstream (rule root_tree), the same
    rooting phylo_mm sees, so it is not re-rooted here.
    """
    tree.ladderize()
    for leaf in tree.get_terminals():
        if leaf.name and leaf.name.startswith("_R_"):
            leaf.name = leaf.name[3:]


def order_display_leaves(leaves: list, pos: dict, alignment: dict, actuation: dict) -> list:
    """Tree order, but runs of identical aligned sequences sorted by % actuation.

    Only identical-sequence (≈zero-length-branch) runs reorder, so the MSA and
    tree geometry do not visibly change while the actuation gradient inside a
    haplotype cluster becomes visible (high at top).
    """
    n = len(leaves)
    order = sorted(range(n), key=lambda i: pos[id(leaves[i])]["y"])
    disp: list[int] = []
    i = 0
    while i < n:
        seq = alignment.get(leaves[order[i]].name)
        j = i + 1
        while seq and j < n and alignment.get(leaves[order[j]].name) == seq:
            j += 1
        run = order[i:j]
        if seq and len(run) > 1:
            run = sorted(
                run,
                key=lambda k: -(actuation.get(leaves[k].name)
                                if actuation.get(leaves[k].name) is not None
                                else float("-inf")),
            )
        disp.extend(run)
        i = j
    return [leaves[k] for k in disp]


def append_panel_rows(
    leaves: list, pos: dict, panel: list[dict] | None,
    alignment: dict, actuation: dict, coverage: dict,
) -> tuple[list, dict, dict, dict, dict[str, str]]:
    """Append synthetic constructs as pseudo-leaf rows at the bottom.

    Constructs carry an aligned sequence (so the MSA shows the edits) but no
    tree branch and no actuation. Each keeps its anchor leaf's x. Returns new
    (leaves, alignment, actuation, coverage, construct_label); the input dicts
    are not modified.
    """
    alignment, actuation, coverage = dict(alignment), dict(actuation), dict(coverage)
    construct_label: dict[str, str] = {}
    if not panel:
        return leaves, alignment, actuation, coverage, construct_label
    augmented: list = list(leaves)
    anchor_x = {lf.name: pos[id(lf)]["x"] for lf in leaves}
    for p in panel:
        obj = _PanelRow(p["name"], p["label"])
        pos[id(obj)] = {"x": anchor_x.get(p["anchor"], 0.0), "y": 0.0}
        alignment[p["name"]] = p["seq"]
        actuation[p["name"]] = None
        coverage[p["name"]] = None
        construct_label[p["name"]] = p["label"]
        augmented.append(obj)
    return augmented, alignment, actuation, coverage, construct_label


def renumber_rows(tree, leaves: list, pos: dict) -> list[float]:
    """Set leaf y to its row index (0..n-1, monotonic) and re-centre internal nodes.

    The Heatmap/Bar traces need a monotonic y, so the display order defines it.
    Returns the per-row y list.
    """
    for row, lf in enumerate(leaves):
        pos[id(lf)]["y"] = float(row)

    def recompute(clade):
        if clade.is_terminal():
            return pos[id(clade)]["y"]
        cys = [recompute(c) for c in clade.clades]
        pos[id(clade)]["y"] = (cys[0] + cys[-1]) / 2.0
        return pos[id(clade)]["y"]

    recompute(tree.root)
    return [pos[id(lf)]["y"] for lf in leaves]


# ── Per-row cell data ─────────────────────────────────────────────────────────

def actuation_cells(names: list[str], acts: list, covs: list) -> tuple[list[str], list[str]]:
    """Per-row actuation bar colour and hover text."""
    colors = [
        ACT_COLORS[min(int(v * 10), 9)] if v is not None else ACT_MISSING
        for v in acts
    ]
    hover = [
        f"<b>{nm}</b><br>Actuation: {f'{v:.1%}' if v is not None else 'missing'}"
        f"<br>Coverage: {c if c is not None else 'missing'}"
        for nm, v, c in zip(names, acts, covs)
    ]
    return colors, hover


def pclai_cells(
    pclai_bed: Path | None, peak_id: str, names: list[str],
) -> tuple[list[str], list[str], bool]:
    """Per-row PCLAI stripe colour + hover, and whether any leaf matched."""
    pclai: dict[str, dict] = {}
    n = len(names)
    if pclai_bed is not None:
        pclai = query_pclai(pclai_bed, peak_id, names)
        n_inherited = sum(1 for v in pclai.values() if v.get("distance_bp", 0) > 0)
        print(
            f"PCLAI ({peak_id}): annotated {len(pclai)}/{n} leaves "
            f"({n_inherited} inherited from nearest ≤ {PCLAI_NEAREST_BP // 1000} kb)",
            file=sys.stderr,
        )
        if not pclai:
            # Usually a coordinate-space or chrom-naming mismatch.
            print(
                f"WARNING: PCLAI ({peak_id}): bed {pclai_bed} supplied but 0/{n} "
                f"leaves matched; the ancestry stripe is blank. Check the peak's "
                f"chrom exists in the BED and leaf names match {{sample}}/h{{hap}}.",
                file=sys.stderr,
            )
    rgbs: list[str] = []
    hover: list[str] = []
    for nm in names:
        info = pclai.get(nm)
        if info is None:
            rgbs.append(PCLAI_MISSING)
            hover.append(f"<b>{nm}</b><br>PCLAI: missing")
        else:
            rgbs.append(f"rgb({info['rgb']})")
            dist = info.get("distance_bp", 0)
            origin = (
                "direct overlap"
                if dist == 0
                else f"<i>nearest window {dist/1000:.1f} kb away</i>"
            )
            hover.append(
                f"<b>{nm}</b><br>"
                f"PCLAI centroid: {info['centroid']}<br>"
                f"segment PC: ({info['pc1']:.3f}, {info['pc2']:.3f})<br>"
                f"confidence: {info['confidence']}<br>"
                f"source: {origin}"
            )
    return rgbs, hover, any(nm in pclai for nm in names)


def msa_window(meta: dict, aln_len: int) -> tuple[bool, int, int]:
    """(re_only, x0, x1): the displayed MSA columns (0-based, inclusive).

    Unless attribution ran on the whole region, the view is clipped to the RE
    columns phylo_mm reports. The heatmap still holds every column; only the
    x-axis range changes, so variant positions keep absolute coordinates.
    """
    lo, hi = meta.get("re_first_col"), meta.get("re_last_col")
    re_only = (
        str(meta.get("attribution_mode", "")).strip() != "region"
        and isinstance(lo, (int, float))
        and isinstance(hi, (int, float))
    )
    if re_only:
        return True, int(lo) - 1, int(hi) - 1
    return False, 0, aln_len - 1


def encode_msa(names: list[str], raw_seqs: list[str], aln_len: int):
    """Encode the MSA for the heatmap.

    Returns (raw_z, consensus, msa_z, msa_cell, msa_hover). raw_z codes A/T/G/C
    as 0-3 and anything else as 5. consensus is the majority non-gap code per
    column (-1 if all gaps). msa_z is 4 where a base matches consensus, 0-3
    where it differs, 5 for gap/N; msa_cell carries the letter at variants only.
    """
    raw_z = [
        [BASE_TO_CODE.get(c.upper(), 5) for c in seq]
        for seq in raw_seqs
    ]

    code_to_char = ["A", "T", "G", "C"]
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

    msa_z: list[list[int]] = []
    msa_cell: list[list[str]] = []
    msa_hover: list[list[str]] = []
    for nm, seq, row_raw in zip(names, raw_seqs, raw_z):
        row_z, row_cell, row_hover = [], [], []
        for col_i, (c, code) in enumerate(zip(seq, row_raw)):
            cu = c.upper()
            cons = consensus[col_i]
            if code == 5:
                row_z.append(5)
                row_cell.append("")
                row_hover.append(f"<b>{nm}</b><br>pos {col_i+1}: {cu} (missing)")
            elif code == cons:
                row_z.append(4)
                row_cell.append("")
                row_hover.append(f"<b>{nm}</b><br>pos {col_i+1}: {cu} (consensus)")
            else:
                row_z.append(code)
                row_cell.append(cu)
                cons_char = code_to_char[cons] if cons >= 0 else "?"
                row_hover.append(
                    f"<b>{nm}</b><br>pos {col_i+1}: {cu} (variant, cons={cons_char})"
                )
        msa_z.append(row_z)
        msa_cell.append(row_cell)
        msa_hover.append(row_hover)
    return raw_z, consensus, msa_z, msa_cell, msa_hover


# ── Layout geometry ───────────────────────────────────────────────────────────

@dataclass
class Layout:
    """Figure size and the paper-coordinate domain of every panel."""
    fig_width: int
    fig_height: int
    label_fs: int
    show_labels: bool
    tree_x: list[float]
    label_x: list[float]
    par_x: list[float]
    pcl_x: list[float]
    act_x: list[float]
    drv_x: list[float]
    msa_x: list[float]
    mark_x: list[float]
    vio_x: list[float]
    logo_y: list[float]
    main_y: list[float]
    vio_y: list[float]


def compute_layout(
    n: int, n_drv: int, msa_disp_n: int, names: list[str],
    has_paralogs: bool, pclai_present: bool,
) -> Layout:
    """Adaptive figure size and panel domains.

    Columns, left to right: TREE | LABEL | PAR | PCL | ACT | DRV | MSA | MARK |
    VIO; PAR/PCL/DRV/MARK/VIO collapse to zero width when their data is absent.
    Height is capped (~2200 px of rows) so big peaks keep a sane aspect ratio;
    width scales with the displayed MSA columns (~4.5 px each, 1500-5200 px).
    Labels are dropped when rows are too short to read, giving the tree the room.
    """
    fig_width = max(1500, min(5200, int(msa_disp_n * 4.5)))
    if n_drv:                       # room for the (narrow, stacked) violin panel
        fig_width = min(6000, fig_width + 300)
    target_h = 2200
    row_h = max(4, min(18, target_h // max(n, 1)))
    label_fs = max(7, min(12, row_h - 1))
    show_labels = label_fs >= 8 and row_h >= 9
    fig_height = max(480, n * row_h + (300 if n_drv else 160))   # +bottom legend

    max_label_chars = max((len(nm) for nm in names), default=10)
    if show_labels:
        label_px = max_label_chars * label_fs * 0.55
        label_frac = min(0.32, max(0.10, label_px / fig_width))
    else:
        label_frac = 0.0
    act_frac = 0.05
    drv_frac = min(0.06, 0.016 * n_drv)   # per-driver carrier columns (left)
    par_frac = 0.018 if has_paralogs else 0.0
    pcl_frac = 0.018 if pclai_present else 0.0
    gap = 0.005
    mark_frac = min(0.02, 0.008 * n_drv)   # one symbol col per driver
    vio_frac = min(0.14, 240 / fig_width) if n_drv else 0.0   # ~2 inch wide
    tree_frac = 0.20
    msa_frac = 1.0 - (label_frac + par_frac + pcl_frac + act_frac + drv_frac
                      + mark_frac + vio_frac + tree_frac + gap * 9)
    cur = 0.0
    tree_x = [cur, cur + tree_frac]
    cur = tree_x[1] + gap
    label_x = [cur, cur + label_frac]
    cur = label_x[1] + gap
    par_x = [cur, cur + par_frac]
    cur = par_x[1] + gap
    pcl_x = [cur, cur + pcl_frac]
    cur = pcl_x[1] + gap
    act_x = [cur, cur + act_frac]
    cur = act_x[1] + gap
    drv_x = [cur, cur + drv_frac]
    cur = drv_x[1] + gap * 2
    if n_drv:
        msa_x = [cur, cur + msa_frac]
        cur = msa_x[1] + gap * 2
        mark_x = [cur, cur + mark_frac]
        cur = mark_x[1] + gap
        vio_x = [cur, 1.00]
    else:
        msa_x = [cur, 1.00]
        mark_x = [0.998, 0.999]   # hidden — no drivers, no shape column/violins
        vio_x = [0.999, 1.00]

    # Vertical: fixed-pixel title strip (~44 px) | logo (~50 px) | main panels.
    title_strip = min(0.10, 44.0 / fig_height)
    logo_band = min(0.08, 50.0 / fig_height)
    logo_y = [1.0 - title_strip - logo_band, 1.0 - title_strip]
    main_y = [0.00, logo_y[0] - 0.005]
    # Violins are ~2-inch boxes stacked from the top; they shrink only on overflow.
    if n_drv:
        plot_h_px = max(400, fig_height - 20 - 170)            # minus t/b margins
        per_vio = 200.0 / plot_h_px
        vio_band = min(main_y[1] - main_y[0], n_drv * per_vio * 1.15)
        vio_y = [main_y[1] - vio_band, main_y[1]]
    else:
        vio_y = main_y
    return Layout(fig_width, fig_height, label_fs, show_labels, tree_x, label_x,
                  par_x, pcl_x, act_x, drv_x, msa_x, mark_x, vio_x, logo_y,
                  main_y, vio_y)


# ── Tree panel ────────────────────────────────────────────────────────────────

def driver_leaf_colors(drivers: list[dict]) -> dict[str, str]:
    """sample_hap → driver colour; a leaf in several drivers takes the smallest."""
    color: dict[str, str] = {}
    size: dict[str, int] = {}
    for d in sorted(drivers, key=lambda d: d["idx"]):
        n_below = d["rec"].get("n_below") or 10**9
        for tip in d["rec"].get("descendants", ()):
            sh = _sample_hap(tip)
            if sh not in size or n_below < size[sh]:
                size[sh] = n_below
                color[sh] = d["color"]
    return color


def add_leader_lines(
    fig, leaves: list, pos: dict, ly: list[float], lay: Layout,
    max_depth: float, leaf_color: dict[str, str],
) -> None:
    """Dotted connectors from each tree tip to its label (or across the figure
    when labels are hidden). Driver carriers get a heavier line in its colour."""
    n = len(leaves)
    x1 = lay.label_x[1] if lay.show_labels else 1.0
    axis_min = -0.01 * max_depth
    axis_max = max_depth * 1.05
    axis_span = axis_max - axis_min
    main_span = lay.main_y[1] - lay.main_y[0]
    for i, leaf in enumerate(leaves):
        tip_x = lay.tree_x[0] + (
            (pos[id(leaf)]["x"] - axis_min) / axis_span
            * (lay.tree_x[1] - lay.tree_x[0])
        )
        # Reversed yaxis: data y=-0.5 is at top, y=n-0.5 is at bottom.
        row_y = lay.main_y[1] - (ly[i] + 0.5) / n * main_span
        c = leaf_color.get(_sample_hap(leaf.name))
        fig.add_shape(
            type="line",
            x0=tip_x, x1=x1,
            y0=row_y, y1=row_y,
            xref="paper", yref="paper",
            line=dict(color=c if c else "#dddddd", width=0.8 if c else 0.4, dash="dot"),
            layer="below",
        )


def add_tree(fig, tree, pos: dict, leaves: list, leaf_text_colors: list[str]) -> None:
    """Grey internal branches, plus terminal branches coloured like the labels
    (one trace per colour)."""
    int_xs, int_ys, term_segs = branch_xy_split(tree.root, pos)
    fig.add_trace(go.Scatter(
        x=int_xs, y=int_ys,
        mode="lines",
        line=dict(color=INTERNAL_BRANCH, width=0.9),
        hoverinfo="skip",
        showlegend=False,
        xaxis="x1", yaxis="y1",
    ))
    terminal_by_color: dict[str, tuple[list, list]] = {}
    for leaf, color in zip(leaves, leaf_text_colors):
        seg = term_segs.get(id(leaf))
        if seg is None:
            continue
        xs, ys = terminal_by_color.setdefault(color, ([], []))
        xs.extend(seg[0] + [None])
        ys.extend(seg[1] + [None])
    for color, (xs, ys) in terminal_by_color.items():
        fig.add_trace(go.Scatter(
            x=xs, y=ys,
            mode="lines",
            line=dict(color=color, width=1.4 if color != LEAF_DEFAULT else 0.9),
            hoverinfo="skip",
            showlegend=False,
            xaxis="x1", yaxis="y1",
        ))


def add_bootstrap(fig, tree, pos: dict) -> None:
    """Bootstrap support labels (≥70) on internal nodes."""
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


class CladeIndex:
    """Parent and tip-set (sample_hap) lookups for every clade, built once."""

    def __init__(self, tree):
        self.tree = tree
        self.tips: dict[int, frozenset] = {}
        for c in tree.find_clades(order="postorder"):
            if c.is_terminal():
                self.tips[id(c)] = frozenset([_sample_hap(c.name)])
            else:
                self.tips[id(c)] = frozenset().union(*(self.tips[id(k)] for k in c.clades))
        parent = {id(k): c for c in tree.find_clades() for k in c.clades}
        # Non-root internal clades in level order, with their parent.
        self.internal = [
            (c, parent[id(c)])
            for c in tree.root.find_clades(order="level")
            if not c.is_terminal() and c is not tree.root
        ]
        self.all_tips = frozenset(_sample_hap(t.name) for t in tree.get_terminals())


def find_branch_endpoints(rec: dict, idx: CladeIndex, pos: dict) -> tuple[float, float, float] | None:
    """Resolve a phylo_mm row to (parent_x, clade_x, clade_y), or None.

    The carrier set can be either half of a bipartition in the drawn tree: take
    the clade whose tips equal it, else the clade equal to its complement and
    use the sibling holding the carriers. Uses the node's own y (the recursive
    rule in compute_positions), not a centroid of the carrier tips.
    """
    want = frozenset(_sample_hap(d) for d in rec.get("descendants", ()))
    if not want:
        return None
    want_complement = idx.all_tips - want
    for clade, parent in idx.internal:
        tips = idx.tips[id(clade)]
        if tips == want:
            return (pos[id(parent)]["x"], pos[id(clade)]["x"], pos[id(clade)]["y"])
        if tips == want_complement:
            for sibling in parent.clades:
                if sibling is not clade and idx.tips[id(sibling)] == want:
                    return (pos[id(parent)]["x"], pos[id(sibling)]["x"], pos[id(sibling)]["y"])
            # Carriers span several siblings (rare on binary trees): use the
            # matched clade with the carriers' centroid y.
            want_ys = [
                pos[id(t)]["y"] for t in idx.tree.get_terminals()
                if _sample_hap(t.name) in want
            ]
            cy = sum(want_ys) / len(want_ys) if want_ys else pos[id(clade)]["y"]
            return (pos[id(parent)]["x"], pos[id(clade)]["x"], cy)
    return None


def add_driver_branches(
    fig, drivers: list[dict], idx: CladeIndex, pos: dict, re_band,
) -> None:
    """Thick coloured branch for each driver whose carriers form a clade.

    The β/p headline is not drawn here (carriers are often not a clade); it sits
    on the driver's MSA column (add_driver_labels), which also owns the legend.
    """
    for gi, drv in enumerate(drivers):
        rep = drv["rec"]
        seg = find_branch_endpoints(rep, idx, pos)
        if seg is None:
            continue
        px, cx, py = seg
        links = drv["links"]
        sorted_sites = (
            [_site_label(s, re_band) for s in sorted(links, key=int)] if links else []
        )
        linked_str_full = ", ".join(sorted_sites) if sorted_sites else "none"
        if len(sorted_sites) <= 4:
            linked_str_short = ", ".join(sorted_sites) if sorted_sites else "no SNP"
        else:
            linked_str_short = (
                ", ".join(sorted_sites[:4]) + f", +{len(sorted_sites) - 4} more"
            )
        beta = rep.get("beta")
        ci_lo, ci_hi = rep.get("ci_lower"), rep.get("ci_upper")
        ci_str = (
            f" [{ci_lo:+.2f}, {ci_hi:+.2f}]"
            if isinstance(ci_lo, float) and isinstance(ci_hi, float)
            else ""
        )
        mi, mo = rep.get("mean_in"), rep.get("mean_out")
        means_str = (
            f"actuation in: {mi:.2f}  vs out: {mo:.2f}<br>"
            if isinstance(mi, float) and isinstance(mo, float)
            else ""
        )
        # Headline = joint (conditional) β and p from the forward selection;
        # the marginal carriers-vs-rest Δ goes alongside in hover.
        head_p = rep.get("p")
        fig.add_trace(go.Scatter(
            x=[px, cx], y=[py, py],
            mode="lines",
            line=dict(color=drv["color"], width=3.0),
            hoverinfo="text",
            hovertext=(
                f"<b>Variant at {linked_str_full}</b>  n_carriers={rep['n_below']}<br>"
                f"{means_str}"
                f"<b>β</b> (joint contribution){ci_str} = {_fmt_beta(beta)}<br>"
                f"<b>p</b> (joint) = {_fmt_p(head_p)}<br>"
                f"<i>marginal Δ (carriers − rest) = {_fmt_beta(_marginal_delta(rep))}, "
                f"BH p_fdr = {_fmt_p(rep.get('p_fdr'))}</i>"
            ),
            name=f"β={_fmt_beta(beta)} p={_fmt_p(head_p)} pos: {linked_str_short}",
            legendgroup=f"sig_group_{gi}",
            showlegend=False,
            xaxis="x1", yaxis="y1",
        ))


# ── Driver panels ─────────────────────────────────────────────────────────────

def add_driver_columns(fig, drivers: list[dict], sh_ly: dict[str, float]) -> None:
    """One column per driver; a filled cell = that haplotype carries it.

    Set-nesting shows only here: an inner driver's filled rows fall inside the
    outer driver's, even when the carriers are not a clade.
    """
    for j, drv in enumerate(drivers):
        cy = []
        for d in drv["rec"].get("descendants", ()):
            sh = _sample_hap(d)
            if sh in sh_ly:
                cy.append(sh_ly[sh])
        if not cy:
            continue
        fig.add_trace(go.Bar(
            x=[0.86] * len(cy), base=[j - 0.43] * len(cy), y=cy,
            orientation="h", width=0.9,
            marker=dict(color=drv["color"], line=dict(width=0)),
            hoverinfo="skip", showlegend=False,
            xaxis="x8", yaxis="y1",
        ))


def add_driver_violins(
    fig, drivers: list[dict], sh_ly: dict[str, float], sh_act: dict[str, float], re_band,
) -> list[str]:
    """Shape column (right of the MSA) + one horizontal split violin per driver.

    A driver's symbol marks its carriers in the shape column and in its violin
    (carriers below the line, comparison group above). The comparison group is
    the local contrast: for a nested driver, the smallest enclosing driver's
    other carriers; otherwise everyone else. This matches the joint β. Returns
    the violin tick labels.
    """
    ticktext: list[str] = []
    all_sh = set(sh_act)
    driver_carriers = [
        {_sample_hap(d) for d in drv["rec"].get("descendants", ())}
        for drv in drivers
    ]
    for j, drv in enumerate(drivers):
        sym, col = drv["symbol"], drv["color"]
        carrier_sh = driver_carriers[j]
        parents = [k for k in range(len(driver_carriers))
                   if k != j and driver_carriers[k] > carrier_sh]
        if parents:
            comp_sh = driver_carriers[
                min(parents, key=lambda k: len(driver_carriers[k]))
            ] - carrier_sh
        else:
            comp_sh = all_sh - carrier_sh
        car_ly = [sh_ly[s] for s in carrier_sh if s in sh_ly]
        rest_ly = [sh_ly[s] for s in comp_sh if s in sh_ly]
        if rest_ly:
            fig.add_trace(go.Scatter(
                x=[j] * len(rest_ly), y=rest_ly, mode="markers",
                marker=dict(symbol="circle-open", color="#cccccc", size=3),
                hoverinfo="skip", showlegend=False, xaxis="x9", yaxis="y1",
            ))
        if car_ly:
            fig.add_trace(go.Scatter(
                x=[j] * len(car_ly), y=car_ly, mode="markers",
                marker=dict(symbol=sym, color=col, size=7,
                            line=dict(width=0.5, color="white")),
                hoverinfo="skip", showlegend=False, xaxis="x9", yaxis="y1",
            ))
        car_x = [sh_act[s] * 100 for s in carrier_sh if s in sh_act]
        rest_x = [sh_act[s] * 100 for s in comp_sh if s in sh_act]
        if rest_x:
            fig.add_trace(go.Violin(
                y=[j] * len(rest_x), x=rest_x, orientation="h",
                side="positive", width=0.9,
                line_color="#bbbbbb", fillcolor="rgba(180,180,180,0.22)",
                points="all", pointpos=0.4, jitter=0.3,
                marker=dict(symbol="circle-open", color="#999999", size=3),
                meanline_visible=True, scalemode="width", spanmode="hard",
                hoverinfo="x", showlegend=False, xaxis="x10", yaxis="y3",
            ))
        if car_x:
            fig.add_trace(go.Violin(
                y=[j] * len(car_x), x=car_x, orientation="h",
                side="negative", width=0.9,
                line_color=col, fillcolor=_hex_rgba(col, 0.25),
                points="all", pointpos=-0.4, jitter=0.3,
                marker=dict(symbol=sym, color=col, size=5),
                meanline_visible=True, scalemode="width", spanmode="hard",
                hoverinfo="x", showlegend=False, xaxis="x10", yaxis="y3",
            ))
        links = drv["links"]
        ticktext.append(_site_label(min(links), re_band) if links else f"d{j}")
        beta = drv["rec"].get("beta")
        bits = []
        # Δ is the local contrast (carriers − comparison), so it agrees with β.
        if car_x and rest_x:
            bits.append(f"Δ={(sum(car_x) / len(car_x) - sum(rest_x) / len(rest_x)):+.0f}%")
        if isinstance(beta, float):
            bits.append(f"β={beta:+.2f}")
        if bits:
            fig.add_annotation(
                x=111, y=j, xref="x10", yref="y3", text="  ".join(bits),
                showarrow=False, font=dict(size=8, color=col), xanchor="right",
            )
    # one-time legend keys: filled symbol = carrier, open ○ = compared-against
    fig.add_trace(go.Scatter(
        x=[None], y=[None], mode="markers",
        marker=dict(symbol="circle", color="#555555", size=7),
        name="violin: driver carriers", legendgroup="vio", showlegend=True,
        legendrank=4000, xaxis="x9", yaxis="y1",
    ))
    fig.add_trace(go.Scatter(
        x=[None], y=[None], mode="markers",
        marker=dict(symbol="circle-open", color="#999999", size=6),
        name="violin: compared against (rest)", legendgroup="vio", showlegend=True,
        legendrank=4001, xaxis="x9", yaxis="y1",
    ))
    return ticktext


# ── Row stripes ───────────────────────────────────────────────────────────────

def add_leaf_labels(fig, names, metas, ly, construct_label, label_fs, colors, hover) -> None:
    """Bare sample-hap labels (construct labels for panel rows), coloured by paralog."""
    display_labels = [
        construct_label.get(nm) or m["sample_hap"]
        for nm, m in zip(names, metas)
    ]
    fig.add_trace(go.Scatter(
        x=[0.0] * len(names), y=ly,
        mode="text",
        text=display_labels,
        textposition="middle right",
        textfont=dict(size=label_fs, color=colors),
        hovertext=hover,
        hoverinfo="text",
        showlegend=False,
        xaxis="x5", yaxis="y1",
    ))


def add_stripe(fig, ly, colors, hover, axis: str) -> None:
    """Full-row colour stripe (paralog / PCLAI).

    Only call when the stripe has data: a zero-width axis with a visible
    full-width Bar makes Plotly paint the whole figure in the bar colour.
    """
    fig.add_trace(go.Bar(
        x=[1] * len(ly), y=ly,
        orientation="h",
        marker=dict(color=colors, line=dict(color="#aaaaaa", width=0.3)),
        width=1.0,
        hovertext=hover,
        hoverinfo="text",
        showlegend=False,
        xaxis=axis, yaxis="y1",
    ))


def add_actuation_bars(fig, ly, colors, hover) -> None:
    fig.add_trace(go.Bar(
        x=[1] * len(ly),
        y=ly,
        orientation="h",
        marker=dict(color=colors, line=dict(color="#aaaaaa", width=0.4)),
        width=0.82,
        hovertext=hover,
        hoverinfo="text",
        showlegend=False,
        xaxis="x2", yaxis="y1",
    ))


def add_legend_keys(fig) -> None:
    """Legend-only traces: actuation bins, FIRE missing, and the two MSA greys."""
    for ri, (color, label) in enumerate(zip(ACT_COLORS, ACT_LABELS)):
        fig.add_trace(go.Bar(
            x=[None], y=[None],
            orientation="h",
            marker=dict(color=color, line=dict(color="#aaaaaa", width=0.4)),
            name=label,
            showlegend=True, legendrank=100 + ri,
        ))
    fig.add_trace(go.Bar(
        x=[None], y=[None],
        orientation="h",
        marker=dict(color=ACT_MISSING, line=dict(color="#aaaaaa", width=0.4)),
        name="FIRE missing",
        showlegend=True, legendrank=120,
    ))
    for color, label in (
        (MSA_INVARIANT, "MSA consensus"),
        (MSA_MISSING, "MSA gap / N"),
    ):
        fig.add_trace(go.Bar(
            x=[None], y=[None],
            orientation="h",
            marker=dict(color=color, line=dict(color="#aaaaaa", width=0.4)),
            name=label,
            legendgroup="msa",
            showlegend=True,
        ))


# ── MSA panel ─────────────────────────────────────────────────────────────────

def add_msa_heatmap(fig, msa_z, ly, msa_cell, msa_hover, aln_len: int) -> None:
    fig.add_trace(go.Heatmap(
        z=msa_z,
        y=ly,
        x=list(range(aln_len)),
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


def add_re_band(fig, re_band: tuple[int, int], y_range: list[float]) -> None:
    """Dashed lines bracketing the RE columns of a flanked alignment."""
    re_lo, re_hi = re_band
    fig.add_trace(go.Scatter(
        x=[re_lo - 0.5, re_lo - 0.5, None, re_hi + 0.5, re_hi + 0.5],
        y=[y_range[1], y_range[0], None, y_range[1], y_range[0]],
        mode="lines",
        line=dict(color="#111111", width=1.2, dash="dash"),
        hoverinfo="text",
        hovertext=f"regulatory element: cols {re_lo+1}–{re_hi+1}",
        name=f"RE (cols {re_lo+1}–{re_hi+1})",
        legendgroup="re",
        showlegend=True,
        xaxis="x3", yaxis="y1",
    ))


def _contiguous_runs(cols: list[int]) -> list[tuple[int, int]]:
    """Sorted columns → [(start, end)] runs of consecutive columns."""
    cols_sorted = sorted(cols)
    runs: list[tuple[int, int]] = []
    run_start = run_end = cols_sorted[0]
    for c0 in cols_sorted[1:]:
        if c0 == run_end + 1:
            run_end = c0
        else:
            runs.append((run_start, run_end))
            run_start = run_end = c0
    runs.append((run_start, run_end))
    return runs


def add_linked_columns(fig, drivers: list[dict], y_range: list[float]) -> None:
    """Outline each driver's linked MSA columns in its colour.

    A column linked to several drivers takes the first (lowest p_fdr) colour.
    Contiguous runs (e.g. a 50-column indel) draw as one box, singletons as one
    dashed line.
    """
    col_to_color: dict[int, str] = {}
    for drv in drivers:
        for s in drv["links"]:
            try:
                c = int(s) - 1
            except (TypeError, ValueError):
                continue
            if c not in col_to_color:
                col_to_color[c] = drv["color"]
    cols_by_color: dict[str, list[int]] = {}
    for c, color in col_to_color.items():
        cols_by_color.setdefault(color, []).append(c)
    for color, cols in cols_by_color.items():
        link_xs: list = []
        link_ys: list = []
        for r_lo, r_hi in _contiguous_runs(cols):
            if r_lo == r_hi:
                link_xs += [r_lo, r_lo, None]
                link_ys += [y_range[1], y_range[0], None]
            else:
                x_lo, x_hi = r_lo - 0.5, r_hi + 0.5
                link_xs += [x_lo, x_hi, x_hi, x_lo, x_lo, None]
                link_ys += [y_range[1], y_range[1], y_range[0], y_range[0], y_range[1], None]
        txt = f"linked variants ({len(cols)}) for sig branch group"
        fig.add_trace(go.Scatter(
            x=link_xs, y=link_ys,
            mode="lines",
            line=dict(color=color, width=1.0, dash="dot"),
            hoverinfo="text",
            hovertext=[txt] * len(link_xs),
            name=f"linked variants ({len(cols)})",
            legendgroup="linked",
            showlegend=True,
            xaxis="x3", yaxis="y1",
        ))


def add_driver_labels(fig, drivers: list[dict], y_range: list[float], re_band) -> None:
    """β/p label above each driver's left-most linked column, plus its legend entry.

    The column always exists, unlike a clade, so it carries the headline. Labels
    are staggered vertically so neighbouring columns do not collide.
    """
    msa_top = y_range[1]   # data-y of the MSA top (reversed axis)
    for i, drv in enumerate(drivers):
        rep = drv["rec"]
        positions = sorted(drv["links"], key=int)
        if not positions:
            continue
        c0 = int(positions[0]) - 1   # 0-based column → x3 coordinate
        beta = _fmt_beta(rep.get("beta"))
        p = _fmt_p(rep.get("p"))
        fig.add_annotation(
            x=c0, y=msa_top,
            xref="x3", yref="y1",
            yanchor="bottom",
            yshift=4 + 14 * (i % 4),
            text=f"<b>β={beta}</b> p={p}",
            showarrow=True, arrowhead=0, arrowwidth=0.8,
            arrowcolor=drv["color"], ax=0, ay=-8,
            font=dict(size=8, color=drv["color"]),
            bgcolor="rgba(255,255,255,0.85)",
            bordercolor=drv["color"], borderwidth=0.5,
        )
        short = _site_label(min(drv["links"]), re_band) if drv["links"] else "driver"
        fig.add_trace(go.Scatter(
            x=[None], y=[None], mode="lines",
            line=dict(color=drv["color"], width=3.0),
            name=f"{short} β={beta}",
            legendgroup=f"sig_group_{i}",
            showlegend=True, legendrank=3000 + i,
            xaxis="x3", yaxis="y1",
        ))


def add_logo(fig, freqs: dict[str, list[float]], aln_len: int) -> None:
    """Stacked per-column base frequencies above the MSA."""
    aln_positions = list(range(aln_len))
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


# ── Layout and header ─────────────────────────────────────────────────────────

def apply_layout(
    fig, lay: Layout, n: int, n_drv: int, max_depth: float,
    msa_x0: int, msa_x1: int, pclai_present: bool, has_paralogs: bool,
    vio_ticktext: list[str],
) -> None:
    """Axes, legend and canvas settings, plus the A/T/G/C key beside the logo."""
    fig.update_layout(
        title=None,   # title is rendered as an annotation so it stays inside the header
        height=lay.fig_height,
        width=lay.fig_width,
        plot_bgcolor="white",
        paper_bgcolor="white",
        bargap=0,
        barmode="overlay",
        # Wipe out the default Plotly template so its heatmap colorscale
        # does not override the explicit colorscale on the MSA trace.
        template="none",
        margin=dict(l=10, r=30, t=20, b=170),
        font=dict(family="Arial"),
        legend=dict(
            title=dict(text="", font=dict(size=10)),
            orientation="h",
            x=0.0, y=-0.085,
            xanchor="left", yanchor="top",
            tracegroupgap=14, entrywidth=0,
            font=dict(size=8),
        ),
        xaxis1=dict(
            domain=lay.tree_x,
            anchor="y1",
            range=[-0.01 * max_depth, max_depth * 1.05],
            title=dict(text="branch length", font=dict(size=10)),
            showgrid=False, zeroline=True, zerolinecolor="#dddddd",
            tickfont=dict(size=9),
        ),
        xaxis5=dict(
            domain=lay.label_x,
            anchor="y1",
            range=[0, 1],
            showticklabels=False, showgrid=False, zeroline=False,
        ),
        xaxis2=dict(
            domain=lay.act_x,
            anchor="y1",
            title=dict(text="% FIRE", font=dict(size=10)),
            showticklabels=False, showgrid=False, zeroline=False,
            range=[0, 1.1],
        ),
        xaxis8=dict(
            domain=lay.drv_x,
            anchor="y1",
            title=dict(text="drivers", font=dict(size=9)) if n_drv else None,
            showticklabels=False, showgrid=False, zeroline=False,
            range=[-0.5, max(0.5, n_drv - 0.5)],
        ),
        xaxis6=dict(
            domain=lay.pcl_x,
            anchor="y1",
            title=dict(text="PCLAI", font=dict(size=9)) if pclai_present else None,
            showticklabels=False, showgrid=False, zeroline=False,
            range=[0, 1.1],
            visible=pclai_present,
        ),
        xaxis7=dict(
            domain=lay.par_x,
            anchor="y1",
            title=dict(text="paralog", font=dict(size=9)) if has_paralogs else None,
            showticklabels=False, showgrid=False, zeroline=False,
            range=[0, 1.1],
            visible=has_paralogs,
        ),
        xaxis3=dict(
            domain=lay.msa_x,
            anchor="y1",
            range=[msa_x0 - 0.5, msa_x1 + 0.5],   # RE-only window unless region mode
            title=dict(text="alignment position", font=dict(size=10)),
            showgrid=False, zeroline=False,
            tickfont=dict(size=9),
        ),
        xaxis4=dict(
            domain=lay.msa_x,
            anchor="y2",
            showticklabels=False,
            showgrid=False, zeroline=False,
            matches="x3",   # zoom/pan stays in sync with MSA x-axis
        ),
        yaxis1=dict(
            domain=lay.main_y,
            range=[n - 0.5, -0.5],   # reversed: leaf 0 at top
            showticklabels=False, showgrid=False, zeroline=False,
        ),
        yaxis2=dict(
            domain=lay.logo_y,
            anchor="x4",
            range=[0, 1.05],
            title=dict(text="freq.", font=dict(size=9)),
            showgrid=False, zeroline=False, tickfont=dict(size=8),
        ),
        # Shape column (carrier markers, leaf-row axis) + violins (% actuation).
        xaxis9=dict(
            domain=lay.mark_x, anchor="y1",
            range=[-0.5, max(0.5, n_drv - 0.5)],
            title=dict(text="in violin", font=dict(size=8)) if n_drv else None,
            showticklabels=False, showgrid=False, zeroline=False,
            visible=bool(n_drv),
        ),
        xaxis10=dict(
            domain=lay.vio_x, anchor="y3", range=[-5, 112],
            title=dict(text="% FIRE", font=dict(size=9)) if n_drv else None,
            tickvals=[0, 50, 100], tickfont=dict(size=8),
            showgrid=True, gridcolor="#eeeeee", zeroline=False,
            visible=bool(n_drv),
        ),
        yaxis3=dict(
            domain=lay.vio_y, anchor="x10",
            range=[max(0.7, n_drv - 0.3), -0.7],   # driver 0 at top, stacked down
            tickvals=list(range(len(vio_ticktext))), ticktext=vio_ticktext,
            tickfont=dict(size=8), showgrid=False, zeroline=False,
            visible=bool(n_drv),
        ),
    )

    for i, nuc in enumerate("ATGC"):
        fig.add_annotation(
            text=f"<b>{nuc}</b>",
            x=1.01,
            y=lay.logo_y[1] - (i + 0.5) * (lay.logo_y[1] - lay.logo_y[0]) / 4,
            xref="paper", yref="paper",
            showarrow=False,
            font=dict(size=10, color=DNA_COLORS[nuc]),
            xanchor="left",
        )


def _h2_summary(meta: dict) -> str | None:
    """Header H²_p text, or None when phylo_mm reported no H² estimate.

    Headline = alpha.V=1 posterior with its 95% HPD; a wide HPD is flagged as
    unidentified, alpha.V=0.25 is shown when it is ≥0.05 lower, and [ESS low]
    marks an unconverged chain. Mode `none` means neither the RE nor the region
    test passed, so no drivers are attributed.
    """
    h2_main = meta.get("h2_alpha1_mean")
    if not isinstance(h2_main, (int, float)):
        return None
    h2_low = meta.get("h2_alpha025_mean")
    h2_lo = meta.get("h2_alpha1_lower")
    h2_hi = meta.get("h2_alpha1_upper")
    conv_flag = " [ESS low]" if str(meta.get("h2_converged", "")).upper() == "FALSE" else ""
    hpd_str, wide_flag = "", ""
    if isinstance(h2_lo, (int, float)) and isinstance(h2_hi, (int, float)):
        hpd_str = f" [{h2_lo:.2f}–{h2_hi:.2f}]"
        if (h2_hi - h2_lo) > 0.5:
            wide_flag = " ⚠ unidentified"
    prior_str = (
        f" (α=0.25→{h2_low:.2f})"
        if isinstance(h2_low, (int, float)) and (h2_main - h2_low) >= 0.05
        else ""
    )
    mode = str(meta.get("attribution_mode", "")).strip()
    if mode == "none":
        def h2_str(prefix):
            m = meta.get(f"{prefix}_mean")
            lo, hi = meta.get(f"{prefix}_lower"), meta.get(f"{prefix}_upper")
            if not isinstance(m, (int, float)):
                return "NA"
            hpd = (f" [{lo:.2f}–{hi:.2f}]"
                   if isinstance(lo, (int, float)) and isinstance(hi, (int, float))
                   else "")
            return f"{m:.2f}{hpd} ⚠"
        return (
            f"<b>H²<sub>p</sub></b> not sequence-driven "
            f"(RE {h2_str('h2_re')}, region {h2_str('h2_region')}){conv_flag}"
        )
    scope = "RE" if mode == "RE" else ("region" if mode == "region" else "")
    lbl = f"H²<sub>p</sub>({scope})" if scope else "H²<sub>p</sub>"
    return f"<b>{lbl}</b>≈{h2_main:.2f}{hpd_str}{prior_str}{wide_flag}{conv_flag}"


def add_header(
    fig, peak_id: str, fig_height: int, *, n: int, n_unique_seqs: int, acts: list,
    tree_length: float, var_sites: int, aln_len: int, re_band, flank: int,
    meta: dict, n_drv: int, paralog_color_map: dict[str, str], has_paralogs: bool,
) -> None:
    """Peak title, the one-line summary under it, and paralog chips (right)."""
    valid_acts = [v for v in acts if v is not None]
    mean_act = sum(valid_acts) / len(valid_acts) if valid_acts else None
    summary_bits = [
        f"<b>n</b>={n}",
        f"<b>unique</b>={n_unique_seqs}",
        f"<b>mean act</b>={mean_act:.1%}" if mean_act is not None else "<b>mean act</b>=NA",
        f"<b>tree length</b>={tree_length:.3f}",
        f"<b>variable sites</b>={var_sites}/{aln_len}",
    ]
    if re_band is not None:
        summary_bits.append(
            f"<b>RE</b>=cols {re_band[0]+1}–{re_band[1]+1} (+{flank}bp flank)"
        )
    h2 = _h2_summary(meta)
    if h2 is not None:
        summary_bits.append(h2)
    summary_bits.append(f"<b>contributing variants</b>={n_drv}")
    fig.add_annotation(
        text=f"<b>{peak_id}</b>",
        x=0.0, y=1.0, xref="paper", yref="paper",
        xanchor="left", yanchor="top",
        showarrow=False,
        font=dict(size=14, family="Arial", color="#111111"),
    )
    fig.add_annotation(
        text="  ·  ".join(summary_bits),
        x=0.0, y=1.0 - (24.0 / fig_height),
        xref="paper", yref="paper",
        xanchor="left", yanchor="top",
        showarrow=False,
        font=dict(size=11, family="Arial", color="#444444"),
    )
    # Paralog chips, right-anchored, 4 per line.
    if has_paralogs:
        items = [
            f"<span style='color:{c}'><b>■</b> {p}</span>"
            for p, c in paralog_color_map.items()
        ]
        rows = ["  ".join(items[i : i + 4]) for i in range(0, len(items), 4)]
        fig.add_annotation(
            text="<br>".join(rows),
            x=1.0, y=1.0,
            xref="paper", yref="paper",
            xanchor="right", yanchor="top",
            showarrow=False,
            font=dict(size=11, family="Arial"),
            align="right",
        )


# ── Main figure builder ───────────────────────────────────────────────────────

def make_figure(
    tree,
    actuation: dict[str, float | None],
    coverage: dict[str, int | None],
    alignment: dict[str, str],
    peak_id: str,
    pclai_bed: Path | None = None,
    phylo_mm: dict | None = None,
    flank: int = 0,
    panel: list[dict] | None = None,
) -> go.Figure:
    """Build the tree | labels | stripes | actuation | MSA | violins figure.

    The tree is ladderized and its ``_R_`` leaf prefixes stripped in place; the
    actuation / coverage / alignment dicts are not modified.
    """
    prepare_tree(tree)
    pos = compute_positions(tree.root)
    leaves = order_display_leaves(list(tree.get_terminals()), pos, alignment, actuation)
    leaves, alignment, actuation, coverage, construct_label = append_panel_rows(
        leaves, pos, panel, alignment, actuation, coverage,
    )
    n = len(leaves)
    ly = renumber_rows(tree, leaves, pos)
    names = [lf.name for lf in leaves]

    # Per-row metadata; constructs get their own label colour and are left out
    # of the unique-sequence count.
    metas = [parse_leaf_meta(nm) for nm in names]
    paralog_color_map, leaf_text_colors, has_paralogs = assign_paralog_colors(metas)
    if construct_label:
        leaf_text_colors = [
            PANEL_LABEL_COLOR if nm in construct_label else c
            for nm, c in zip(names, leaf_text_colors)
        ]
    n_unique_seqs = len({
        alignment.get(nm, "")
        for nm in names
        if alignment.get(nm) and nm not in construct_label
    })
    acts = [actuation.get(nm) for nm in names]
    act_colors, act_hover = actuation_cells(names, acts, [coverage.get(nm) for nm in names])
    pclai_rgbs, pclai_hover, pclai_present = pclai_cells(pclai_bed, peak_id, names)

    aln_len = max((len(s) for s in alignment.values()), default=0)
    raw_seqs = [alignment.get(nm, "-" * aln_len) for nm in names]
    meta = (phylo_mm or {}).get("meta", {})
    re_only_msa, msa_x0, msa_x1 = msa_window(meta, aln_len)
    drivers = driver_groups(phylo_mm)
    n_drv = len(drivers)
    re_band = re_band_columns(alignment, names, flank, aln_len)
    raw_z, consensus, msa_z, msa_cell, msa_hover = encode_msa(names, raw_seqs, aln_len)

    lay = compute_layout(n, n_drv, max(1, msa_x1 - msa_x0 + 1), names,
                         has_paralogs, pclai_present)
    max_depth = max((pos[id(lf)]["x"] for lf in leaves), default=1.0)
    y_range = [n - 0.5, -0.5]   # reversed: leaf 0 at top

    # Trace order is fixed: later traces draw on top.
    fig = go.Figure()
    add_leader_lines(fig, leaves, pos, ly, lay, max_depth, driver_leaf_colors(drivers))
    add_tree(fig, tree, pos, leaves, leaf_text_colors)
    vio_ticktext: list[str] = []
    if drivers:
        sh_ly = {
            _sample_hap(names[i]): ly[i]
            for i in range(n) if names[i] not in construct_label
        }
        sh_act = {
            _sample_hap(names[i]): acts[i]
            for i in range(n)
            if names[i] not in construct_label and acts[i] is not None
        }
        add_driver_columns(fig, drivers, sh_ly)
        vio_ticktext = add_driver_violins(fig, drivers, sh_ly, sh_act, re_band)
    add_driver_branches(fig, drivers, CladeIndex(tree), pos, re_band)
    add_bootstrap(fig, tree, pos)
    if lay.show_labels:
        add_leaf_labels(fig, names, metas, ly, construct_label, lay.label_fs,
                        leaf_text_colors, act_hover)
    if has_paralogs:
        paralog_hover = [
            f"<b>{nm}</b><br>paralog: {m['paralog'] or '—'}"
            for nm, m in zip(names, metas)
        ]
        add_stripe(fig, ly, leaf_text_colors, paralog_hover, "x7")
    if pclai_present:
        add_stripe(fig, ly, pclai_rgbs, pclai_hover, "x6")
    add_actuation_bars(fig, ly, act_colors, act_hover)
    add_legend_keys(fig)
    add_msa_heatmap(fig, msa_z, ly, msa_cell, msa_hover, aln_len)
    # RE brackets only when the view is not already clipped to the RE.
    if re_band is not None and not re_only_msa:
        add_re_band(fig, re_band, y_range)
    add_linked_columns(fig, drivers, y_range)
    add_driver_labels(fig, drivers, y_range, re_band)
    add_logo(fig, logo_freqs(raw_z, aln_len), aln_len)
    apply_layout(fig, lay, n, n_drv, max_depth, msa_x0, msa_x1,
                 pclai_present, has_paralogs, vio_ticktext)
    var_sites = sum(1 for col in range(aln_len) if any(
        row[col] != consensus[col] and row[col] < 4 for row in raw_z
    ))
    add_header(
        fig, peak_id, lay.fig_height, n=n, n_unique_seqs=n_unique_seqs, acts=acts,
        tree_length=max_depth, var_sites=var_sites, aln_len=aln_len,
        re_band=re_band, flank=flank, meta=meta, n_drv=n_drv,
        paralog_color_map=paralog_color_map, has_paralogs=has_paralogs,
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
    parser.add_argument("--html",      required=False, default=None, type=Path,
                        help="Output HTML path")
    parser.add_argument("--pdf",       required=True, type=Path,
                        help="Output PDF path (via kaleido)")
    parser.add_argument("--peak-id",   type=str, default=None)
    parser.add_argument(
        "--flank", type=int, default=0,
        help="bp of flank that add_flanks added each side of the RE. When >0 "
             "the RE sub-span is bracketed in the MSA and noted in the header.",
    )
    parser.add_argument(
        "--pclai-bed", type=Path, default=None,
        help="Optional bgzipped + tabix-indexed HPRC PCLAI CHM13-coord BED. "
             "When supplied, draws a per-leaf ancestry stripe colored by HPRC's "
             "pre-baked PCLAI palette.",
    )
    parser.add_argument(
        "--phylo-mm", type=Path, default=None,
        help="Optional per-peak phylo_mm TSV (workflow/scripts/phylo_mm.R). "
             "Adds H²_p and the driver count to the header summary; per-variant "
             "drivers (forward-selection WLS, is_split) are annotated on the "
             "tree with their MSA columns highlighted.",
    )
    parser.add_argument(
        "--panel", type=Path, default=None,
        help="Optional make_mutagenesis_panel --aligned-out FASTA. Its synthetic "
             "constructs render as grey pseudo-leaf rows (MSA + edits, no branch, "
             "no actuation) next to the carrier leaf each was built from, for "
             "visually validating the single-change panel.",
    )
    args = parser.parse_args()

    tree                = Phylo.read(str(args.tree), "newick")
    actuation, coverage = load_actuation(args.actuation)
    alignment           = load_alignment(args.alignment)
    peak_id             = args.peak_id or args.tree.stem

    phylo_mm = load_phylo_mm(args.phylo_mm) if args.phylo_mm else {"meta": {}, "branches": []}
    panel = load_panel(args.panel) if args.panel else None

    fig = make_figure(
        tree, actuation, coverage, alignment, peak_id,
        pclai_bed=args.pclai_bed,
        phylo_mm=phylo_mm,
        flank=args.flank,
        panel=panel,
    )

    # HTML is optional (interactive but ~4 MB each) — only when --html is given.
    if args.html is not None:
        args.html.parent.mkdir(parents=True, exist_ok=True)
        fig.write_html(
            str(args.html),
            include_plotlyjs="cdn",
            config={
                "toImageButtonOptions": {"format": "svg", "filename": peak_id},
                "displayModeBar": True,
            },
        )
        print(f"Saved {args.html}", file=sys.stderr)

    # PDF via kaleido. Fail loudly: snakemake declares `pdf` as a required
    # output, so a silent PDF failure produces a confusing "missing outputs"
    # error instead of the actual kaleido stack trace.
    fig.write_image(str(args.pdf), format="pdf", scale=2)
    print(f"Saved {args.pdf}", file=sys.stderr)


if __name__ == "__main__":
    main()
