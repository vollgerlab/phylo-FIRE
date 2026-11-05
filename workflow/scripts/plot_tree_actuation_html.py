#!/usr/bin/env python3
"""Interactive phylogenetic tree + FIRE actuation + MSA visualization via Plotly.

Layout (left to right):
  tree | leaf labels | PCLAI ancestry stripe | FIRE % actuation | MSA
                                              with sequence logo above

Always writes HTML. Also writes PDF via kaleido (already in pixi deps).

The PCLAI stripe uses HPRC R2's pre-baked PCLAI palette (col 9 of the bgzipped
+ tabix-indexed CHM13-coord BED at resources/pclai/hprc_r2.bed.gz). The
peak_id is parsed as ``{chrom}_{start}_{end}`` and a single tabix query
returns one row per sample-haplotype overlapping the peak; each row's sample/hap
identity is taken from BED col 4. Leaves use ``{sample}_{hap}`` directly
(union mode) or ``{locus}__{paralog}__{sample}_{hap}`` (paralog mode); only
the trailing ``__``-segment is matched against PCLAI.

Peaks whose peak_id chrom isn't a CHM13 chromosome (e.g.
``chr8_MATERNAL_8587817_8588089`` HG002-only consensus) silently fall through
to no PCLAI stripe. The stripe is also omitted when --pclai-bed is not supplied.

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
from pathlib import Path

from Bio import Phylo
import plotly.graph_objects as go

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


# ── PCLAI lookup ──────────────────────────────────────────────────────────────

PCLAI_MISSING = "#d0d0d0"   # same grey used for missing actuation
PCLAI_NEAREST_BP = 1_000_000  # max gap to inherit ancestry from a neighbouring window


def _leaf_to_sample_hap(leaf_name: str) -> str:
    """Extract ``{sample}_{hap}`` from a leaf name.

    Union-mode leaves are already in that form (``HG02886_1``).
    Paralog-mode leaves carry ``{locus}__{paralog}__{sample}_{hap}``;
    we keep only the trailing ``__``-segment.
    """
    return leaf_name.split("__")[-1]


def _parse_peak_coords(s: str) -> tuple[str, int, int] | None:
    """Parse ``{chrom}_{start}_{end}`` (e.g. ``chr1_120083232_120083371``)."""
    parts = s.rsplit("_", 2)
    if len(parts) != 3:
        return None
    try:
        return parts[0], int(parts[1]), int(parts[2])
    except ValueError:
        return None


def _row_to_info(f: list[str], distance_bp: int = 0) -> dict:
    """Convert a parsed tabix row into the per-leaf info dict.

    ``distance_bp`` is 0 for rows that directly overlap the peak, and >0 when
    the row was picked by the nearest-window fallback.
    """
    try:
        pc_str = f[3].rsplit("_", 1)[-1].strip("()")
        pc1, pc2 = (float(x) for x in pc_str.split(","))
    except (IndexError, ValueError):
        pc1, pc2 = float("nan"), float("nan")
    return {
        "rgb": f[8],
        "centroid": f[9],
        "pc1": pc1,
        "pc2": pc2,
        "confidence": int(f[4]),
        "distance_bp": distance_bp,
    }


def _overlap_score(rs: int, re_: int, start: int, end: int) -> tuple[int, int] | None:
    """(cost, distance_bp) for exact overlap: larger overlap wins; no overlap drops."""
    overlap = max(0, min(re_, end) - max(rs, start))
    return None if overlap <= 0 else (-overlap, 0)


def _nearest_score(rs: int, re_: int, start: int, end: int) -> tuple[int, int] | None:
    """(cost, distance_bp) for the nearest window within PCLAI_NEAREST_BP (0 = overlap)."""
    if re_ <= start:
        dist = start - re_
    elif rs >= end:
        dist = rs - end
    else:
        dist = 0
    return None if dist > PCLAI_NEAREST_BP else (dist, dist)


def _match_pclai(rows, start, end, leaf_names, score) -> dict[str, dict]:
    """Pick the lowest-cost tabix row per sample-hap and fan it out to its leaves.

    ``score(rs, re_, start, end)`` returns ``(cost, distance_bp)`` or None to
    drop the row; ties keep the first row seen.
    """
    leaves_by_sample_hap: dict[str, list[str]] = {}
    for nm in leaf_names:
        leaves_by_sample_hap.setdefault(_leaf_to_sample_hap(nm), []).append(nm)

    best: dict[str, tuple[int, int, list[str]]] = {}
    for raw in rows:
        f = raw.split("\t")
        try:
            sample, hap_str = f[3].split("/")[:2]
        except (IndexError, ValueError):
            continue
        sample_hap = f"{sample}_{hap_str.lstrip('h')}"
        if sample_hap not in leaves_by_sample_hap:
            continue
        sc = score(int(f[1]), int(f[2]), start, end)
        if sc is None:
            continue
        prev = best.get(sample_hap)
        if prev is None or sc[0] < prev[0]:
            best[sample_hap] = (sc[0], sc[1], f)

    out: dict[str, dict] = {}
    for sample_hap, (_, dist, f) in best.items():
        info = _row_to_info(f, distance_bp=dist)
        for leaf in leaves_by_sample_hap[sample_hap]:
            out[leaf] = info
    return out


def _fetch(tabix, chrom: str, start: int, end: int) -> list[str]:
    """Tabix rows in [start, end), or [] when the chrom is absent."""
    try:
        return list(tabix.fetch(chrom, start, end))
    except (ValueError, OSError):
        return []


def query_pclai(pclai_bed: Path, peak_id: str, leaf_names: list[str]) -> dict[str, dict]:
    """Per-leaf PCLAI lookup in CHM13 space.

    Union peaks parse ``peak_id`` as ``{chrom}_{start}_{end}``; paralog groups
    take each leaf's sub-peak from the first ``__``-segment of its name. Each
    leaf gets the window with the largest overlap, else the nearest window
    within PCLAI_NEAREST_BP (``distance_bp`` > 0).

    Returns ``{leaf_name: {"rgb", "centroid", "pc1", "pc2", "confidence",
    "distance_bp"}}``; unannotated leaves are absent.
    """
    import pysam   # lazy: only required when PCLAI is wired in

    parsed = _parse_peak_coords(peak_id)
    if parsed is not None:
        intervals = [(*parsed, leaf_names)]
    else:
        by_peak: dict[tuple[str, int, int], list[str]] = {}
        for nm in leaf_names:
            coords = _parse_peak_coords(nm.split("__", 1)[0])
            if coords:
                by_peak.setdefault(coords, []).append(nm)
        intervals = [(*k, v) for k, v in by_peak.items()]

    out: dict[str, dict] = {}
    tabix = pysam.TabixFile(str(pclai_bed))
    try:
        for chrom, start, end, leaves in intervals:
            rows = _fetch(tabix, chrom, start, end)
            if rows:
                out.update(_match_pclai(rows, start, end, leaves, _overlap_score))
        for chrom, start, end, leaves in intervals:
            unannot = [nm for nm in leaves if nm not in out]
            if not unannot:
                continue
            rows = _fetch(tabix, chrom, max(0, start - PCLAI_NEAREST_BP), end + PCLAI_NEAREST_BP)
            if rows:
                out.update(_match_pclai(rows, start, end, unannot, _nearest_score))
    finally:
        tabix.close()
    return out


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
        try: rec["branch_id"] = int(rec["branch_id"])
        except (KeyError, ValueError): pass
        try: rec["n_below"]   = int(rec["n_below"])
        except (KeyError, ValueError): pass
        for k in ("mean_in", "mean_out", "t_stat", "p_raw", "p_fdr",
                  "lasso_coef", "beta", "se",
                  "ci_lower", "ci_upper", "p"):
            raw_val = rec.get(k, "")
            if raw_val in ("", "NA"):
                rec[k] = None
            else:
                try: rec[k] = float(raw_val)
                except ValueError: rec[k] = None
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


# ── Main figure builder ───────────────────────────────────────────────────────

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
    # The tree is already midpoint-rooted upstream (rule root_tree) — the SAME
    # rooting phylo_mm sees — so we do NOT re-root here. Just ladderize for a
    # tidy display order (ladderize doesn't change the rooting or clades).
    tree.ladderize()

    # MAFFT --adjustdirection prefixes reversed-complemented leaf names with "_R_".
    # IQ-TREE inherits those names; strip the prefix so they match the actuation TSV.
    for _leaf in tree.get_terminals():
        if _leaf.name and _leaf.name.startswith("_R_"):
            _leaf.name = _leaf.name[3:]

    pos    = compute_positions(tree.root)
    leaves = list(tree.get_terminals())
    n      = len(leaves)

    # Display order: tree order (by leaf y), but within each contiguous block of
    # IDENTICAL aligned sequences, order rows by % actuation (high at top) so the
    # actuation gradient inside a haplotype cluster is visible. We reorder the
    # leaves list itself and renumber y to 0..n-1 — kept MONOTONIC because the
    # MSA / actuation go.Heatmap traces (y=ly) render gaps for a non-monotonic y.
    # Only identical-sequence (≈zero-length-branch) runs reorder, so the MSA and
    # tree geometry don't visibly change. Internal-node y is then recomputed.
    _ord = sorted(range(n), key=lambda i: pos[id(leaves[i])]["y"])
    _disp: list[int] = []
    _i = 0
    while _i < n:
        _seq = alignment.get(leaves[_ord[_i]].name)
        _j = _i + 1
        while _seq and _j < n and alignment.get(leaves[_ord[_j]].name) == _seq:
            _j += 1
        run = _ord[_i:_j]
        if _seq and len(run) > 1:
            run = sorted(
                run,
                key=lambda k: -(actuation.get(leaves[k].name)
                                if actuation.get(leaves[k].name) is not None
                                else float("-inf")),
            )
        _disp.extend(run)
        _i = _j
    leaves = [leaves[k] for k in _disp]

    # Optional: append synthetic constructs (mutagenesis / plasmid) as pseudo-leaf
    # rows in a BLOCK at the BOTTOM of the MSA. They carry an aligned sequence (so
    # the MSA shows the edits) but NO tree branch and NO actuation — grey rows for
    # validating the panel. They keep their anchor leaf's x (an unused tree pos);
    # tree branches read the SAME pos[]["y"], so the tree stays in sync and just
    # ends above the construct block.
    construct_names: set[str] = set()
    construct_label: dict[str, str] = {}
    if panel:
        augmented: list = list(leaves)          # every real leaf first
        anchor_x = {lf.name: pos[id(lf)]["x"] for lf in leaves}
        for p in panel:
            obj = _PanelRow(p["name"], p["label"])
            pos[id(obj)] = {"x": anchor_x.get(p["anchor"], 0.0), "y": 0.0}
            alignment[p["name"]] = p["seq"]
            actuation[p["name"]] = None
            coverage[p["name"]] = None
            construct_names.add(p["name"])
            construct_label[p["name"]] = p["label"]
            augmented.append(obj)               # constructs grouped at the bottom
        leaves = augmented
    n = len(leaves)

    for _row, _lf in enumerate(leaves):
        pos[id(_lf)]["y"] = float(_row)

    def _recompute_y(clade):
        if clade.is_terminal():
            return pos[id(clade)]["y"]
        cys = [_recompute_y(c) for c in clade.clades]
        pos[id(clade)]["y"] = (cys[0] + cys[-1]) / 2.0
        return pos[id(clade)]["y"]
    _recompute_y(tree.root)

    names  = [lf.name for lf in leaves]
    ly     = [pos[id(lf)]["y"] for lf in leaves]   # 0 … n-1, monotonic

    # Per-leaf metadata + colour assignment (paralog-aware).
    metas = [parse_leaf_meta(nm) for nm in names]
    paralog_color_map, leaf_text_colors, has_paralogs = assign_paralog_colors(metas)
    if construct_names:
        # constructs get a distinct label colour and are excluded from the
        # real-haplotype unique-sequence stat.
        leaf_text_colors = [
            PANEL_LABEL_COLOR if nm in construct_names else c
            for nm, c in zip(names, leaf_text_colors)
        ]
    n_unique_seqs = len({
        alignment.get(nm, "")
        for nm in names
        if alignment.get(nm) and nm not in construct_names
    })

    # Actuation
    acts   = [actuation.get(nm) for nm in names]
    covs   = [coverage.get(nm) for nm in names]
    colors = [
        ACT_COLORS[min(int(v * 10), 9)] if v is not None else ACT_MISSING
        for v in acts
    ]
    act_hover = [
        f"<b>{nm}</b><br>Actuation: {f'{v:.1%}' if v is not None else 'missing'}"
        f"<br>Coverage: {c if c is not None else 'missing'}"
        for nm, v, c in zip(names, acts, covs)
    ]

    # PCLAI ancestry stripe (only when --pclai-bed was supplied AND at least
    # one leaf has a hit). Uses HPRC's pre-baked itemRgb so colors match
    # the official PCLAI palette. Query happens here, after _R_-prefix
    # stripping, so leaf names match the {sample}/h{hap} keys in BED col 4.
    pclai: dict[str, dict] = {}
    if pclai_bed is not None:
        pclai = query_pclai(pclai_bed, peak_id, names)
        n_inherited = sum(1 for v in pclai.values() if v.get("distance_bp", 0) > 0)
        print(
            f"PCLAI ({peak_id}): annotated {len(pclai)}/{n} leaves "
            f"({n_inherited} inherited from nearest ≤ {PCLAI_NEAREST_BP // 1000} kb)",
            file=sys.stderr,
        )
        if not pclai:
            # The BED exists but no leaf matched — usually a coordinate-space or
            # chrom-naming mismatch rather than a genuine absence of ancestry.
            print(
                f"WARNING: PCLAI ({peak_id}): bed {pclai_bed} supplied but 0/{n} "
                f"leaves matched; the ancestry stripe is blank. Check the peak's "
                f"chrom exists in the BED and leaf names match {{sample}}/h{{hap}}.",
                file=sys.stderr,
            )
    pclai_rgbs: list[str] = []
    pclai_hover: list[str] = []
    for nm in names:
        info = pclai.get(nm)
        if info is None:
            pclai_rgbs.append(PCLAI_MISSING)
            pclai_hover.append(f"<b>{nm}</b><br>PCLAI: missing")
        else:
            pclai_rgbs.append(f"rgb({info['rgb']})")
            dist = info.get("distance_bp", 0)
            origin = (
                "direct overlap"
                if dist == 0
                else f"<i>nearest window {dist/1000:.1f} kb away</i>"
            )
            pclai_hover.append(
                f"<b>{nm}</b><br>"
                f"PCLAI centroid: {info['centroid']}<br>"
                f"segment PC: ({info['pc1']:.3f}, {info['pc2']:.3f})<br>"
                f"confidence: {info['confidence']}<br>"
                f"source: {origin}"
            )
    pclai_present = any(nm in pclai for nm in names)

    # MSA — three-way encoding: variable ATGC | invariant | missing
    aln_len = max((len(s) for s in alignment.values()), default=0)
    raw_seqs = [alignment.get(nm, "-" * aln_len) for nm in names]

    # MSA display window. The flanks only ever matter for (a) building the tree
    # and (b) region-mode drivers that sit outside the RE. So unless attribution
    # expanded to the whole region, clip the displayed MSA to the RE band — the
    # tree is already built, and a 950-column flanked alignment collapses to a
    # readable ~400-column view. The heatmap still holds every column; we just
    # set the x-axis range (so variant positions keep their absolute coords).
    _meta0 = (phylo_mm or {}).get("meta", {})
    _re_lo_m, _re_hi_m = _meta0.get("re_first_col"), _meta0.get("re_last_col")
    re_only_msa = (
        str(_meta0.get("attribution_mode", "")).strip() != "region"
        and isinstance(_re_lo_m, (int, float))
        and isinstance(_re_hi_m, (int, float))
    )
    if re_only_msa:
        msa_x0, msa_x1 = int(_re_lo_m) - 1, int(_re_hi_m) - 1   # 0-based inclusive
    else:
        msa_x0, msa_x1 = 0, aln_len - 1
    msa_disp_n = max(1, msa_x1 - msa_x0 + 1)

    # Contributing variants ("drivers") — counted early so the layout can
    # reserve a per-driver carrier column (filled = the haplotype carries that
    # driver's allele). This is the ONLY place set-nesting is visible: when an
    # inner driver's carriers are a subset of an outer driver's, the inner
    # column's filled cells line up inside the outer's — which the tree can't
    # show because the carriers needn't be a clade.
    drivers = driver_groups(phylo_mm)
    n_drv = len(drivers)

    # Regulatory-element band (only when the alignment carries flanks).
    re_band = re_band_columns(alignment, names, flank, aln_len)

    def _site_label(s) -> str:
        """Label a linked site by RE-relative position when flanks are present.

        phylo_mm reports linked sites as 1-based alignment columns; subtracting
        the RE's 0-based start column yields the 1-based position within the RE
        (so `RE13` reads as "13 bp into the regulatory element", not an opaque
        flanked-alignment column). Without flanks the column already is the RE
        position, so it's shown as-is.
        """
        if re_band is None:
            return str(s)
        try:
            return f"RE{int(s) - re_band[0]}"
        except (TypeError, ValueError):
            return str(s)

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

    # ── Adaptive figure dimensions ────────────────────────────────────────────
    # Cap total figure height so the aspect ratio stays sane for large peaks
    # (paralog groups can have hundreds of leaves; without a cap a 276-leaf peak
    # produces a 3500 px tall figure, which fit-to-width viewers compress to
    # the point that 0.9 px tree lines are sub-pixel and disappear).
    # Scale width with the alignment length so each MSA column keeps enough
    # horizontal room to read its variant marks (the MSA panel is ~55% of the
    # figure, so ~4.5 px/column total → ~2.5 px in the MSA). Floor at 1500 for
    # short RE-only peaks, cap so very long flanked/paralog alignments don't
    # produce an absurd canvas.
    fig_width   = max(1500, min(5200, int(msa_disp_n * 4.5)))
    if n_drv:                       # room for the (narrow, stacked) violin panel
        fig_width = min(6000, fig_width + 300)
    target_h    = 2200
    row_h       = max(4, min(18, target_h // max(n, 1)))
    label_fs    = max(7, min(12, row_h - 1))
    show_labels = label_fs >= 8 and row_h >= 9
    fig_height  = max(480, n * row_h + (300 if n_drv else 160))   # +bottom legend

    # ── Adaptive horizontal layout ────────────────────────────────────────────
    # Horizontal columns:  tree | leaf labels | PCLAI stripe | actuation | MSA
    # When the row is too short to fit readable label text, the label column
    # collapses to zero so the tree gets the reclaimed width; identity is
    # available via hover on the PCLAI / actuation / MSA columns.
    max_label_chars = max((len(nm) for nm in names), default=10)
    if show_labels:
        label_px   = max_label_chars * label_fs * 0.55
        label_frac = min(0.32, max(0.10, label_px / fig_width))
    else:
        label_px   = 0.0
        label_frac = 0.0
    act_frac        = 0.05
    drv_frac        = min(0.06, 0.016 * n_drv)   # per-driver carrier columns (left)
    par_frac        = 0.018 if has_paralogs else 0.0
    pcl_frac        = 0.018 if pclai_present else 0.0
    gap             = 0.005
    # Right-of-MSA driver panels: a per-driver SHAPE column (which haps feed each
    # violin) + the violins (carrier vs rest % actuation). Zero-width when there
    # are no drivers.
    mark_frac       = min(0.02, 0.008 * n_drv)   # tight: one symbol col per driver
    vio_frac        = min(0.14, 240 / fig_width) if n_drv else 0.0   # ~2 inch wide
    tree_frac       = 0.20
    msa_frac        = 1.0 - (label_frac + par_frac + pcl_frac + act_frac + drv_frac
                             + mark_frac + vio_frac + tree_frac + gap * 9)
    # Slide order: TREE | LABEL | PAR | PCL | ACT | DRV | MSA | MARK | VIO.
    # PAR/PCL/DRV/MARK/VIO collapse to zero width when their data isn't present;
    # traces stay in the figure to keep trace indices stable.
    cur = 0.0
    TREE_X  = [cur, cur + tree_frac];   cur = TREE_X[1]  + gap
    LABEL_X = [cur, cur + label_frac];  cur = LABEL_X[1] + gap
    PAR_X   = [cur, cur + par_frac];    cur = PAR_X[1]   + gap
    PCL_X   = [cur, cur + pcl_frac];    cur = PCL_X[1]   + gap
    ACT_X   = [cur, cur + act_frac];    cur = ACT_X[1]   + gap
    DRV_X   = [cur, cur + drv_frac];    cur = DRV_X[1]   + gap * 2
    if n_drv:
        MSA_X   = [cur, cur + msa_frac]; cur = MSA_X[1]  + gap * 2
        MARK_X  = [cur, cur + mark_frac]; cur = MARK_X[1] + gap
        VIO_X   = [cur, 1.00]
    else:
        MSA_X   = [cur, 1.00]
        MARK_X  = [0.998, 0.999]   # hidden — no drivers, no shape column/violins
        VIO_X   = [0.999, 1.00]

    # Vertical layout: title strip (top, ~40 px) | freq chart (~50 px) | main.
    # The title strip stays a fixed pixel height regardless of figure tallness,
    # so the title + paralog chips never collide with the freq chart bars.
    TITLE_STRIP_FRAC = min(0.10, 44.0 / fig_height)
    LOGO_BAND_FRAC   = min(0.08, 50.0 / fig_height)
    LOGO_Y           = [1.0 - TITLE_STRIP_FRAC - LOGO_BAND_FRAC, 1.0 - TITLE_STRIP_FRAC]
    MAIN_Y           = [0.00, LOGO_Y[0] - 0.005]
    # Violins are fixed ~2-inch boxes stacked from the top of the right column;
    # they only shrink when n_drv would overflow the available height.
    if n_drv:
        _plot_h_px = max(400, fig_height - 20 - 170)            # minus t/b margins
        _per_vio   = 200.0 / _plot_h_px                          # ~2 inch each
        _vio_band  = min(MAIN_Y[1] - MAIN_Y[0], n_drv * _per_vio * 1.15)
        VIO_Y3     = [MAIN_Y[1] - _vio_band, MAIN_Y[1]]
    else:
        VIO_Y3     = MAIN_Y

    max_depth = max((pos[id(lf)]["x"] for lf in leaves), default=1.0)

    # Shared y-range — explicit to eliminate auto-padding that breaks alignment
    Y_RANGE = [n - 0.5, -0.5]   # reversed: leaf 0 at top

    fig = go.Figure()

    # ── 0a. Leaf → driver colour for leader lines ────────────────────────────
    # A leaf in several drivers takes the smallest one's colour.
    leaf_clade_color: dict[str, str] = {}
    leaf_clade_size: dict[str, int] = {}
    for _d in sorted(drivers, key=lambda d: d["idx"]):
        _color = _d["color"]
        _n = _d["rec"].get("n_below") or 10**9
        for _tip in _d["rec"].get("descendants", ()):
            _sh = parse_leaf_meta(_tip)["sample_hap"]
            if _sh not in leaf_clade_size or _n < leaf_clade_size[_sh]:
                leaf_clade_size[_sh] = _n
                leaf_clade_color[_sh] = _color

    # ── 0b. Leaf-to-label leader lines ───────────────────────────────────────
    # Light dotted horizontal connectors from each tree tip to the sample-name
    # column, so the eye can trace each tip to its label. They stop at the names
    # (LABEL_X[1]) rather than running on through the PCLAI / actuation / MSA
    # panels. When labels are hidden (tiny rows) there are no names to land on,
    # so the lines extend the full width as the only row guide.
    # Leaves carrying a driver get a heavier line in its colour; others stay gray.
    _leader_x1 = LABEL_X[1] if show_labels else 1.0
    _tree_axis_min = -0.01 * max_depth
    _tree_axis_max = max_depth * 1.05
    _tree_axis_span = _tree_axis_max - _tree_axis_min
    _main_span = MAIN_Y[1] - MAIN_Y[0]
    for _i, _leaf in enumerate(leaves):
        _tip_axis_x = pos[id(_leaf)]["x"]
        _tip_paper_x = TREE_X[0] + (
            (_tip_axis_x - _tree_axis_min) / _tree_axis_span
            * (TREE_X[1] - TREE_X[0])
        )
        # Reversed yaxis: data y=-0.5 is at top, y=n-0.5 is at bottom.
        _row_paper_y = MAIN_Y[1] - (ly[_i] + 0.5) / n * _main_span
        _sample_hap = parse_leaf_meta(_leaf.name)["sample_hap"]
        _clade_color = leaf_clade_color.get(_sample_hap)
        _line_color = _clade_color if _clade_color else "#dddddd"
        _line_width = 0.8 if _clade_color else 0.4
        fig.add_shape(
            type="line",
            x0=_tip_paper_x, x1=_leader_x1,
            y0=_row_paper_y, y1=_row_paper_y,
            xref="paper", yref="paper",
            line=dict(color=_line_color, width=_line_width, dash="dot"),
            layer="below",
        )

    # ── 1. Tree branches ──────────────────────────────────────────────────────
    # Internal segments stay grey; terminal segments are coloured by paralog so
    # the tree carries the same colour signal as the leaf labels.
    int_xs, int_ys, term_segs = branch_xy_split(tree.root, pos)
    fig.add_trace(go.Scatter(
        x=int_xs, y=int_ys,
        mode="lines",
        line=dict(color=INTERNAL_BRANCH, width=0.9),
        hoverinfo="skip",
        showlegend=False,
        xaxis="x1", yaxis="y1",
    ))
    # Group terminal branches by paralog so we emit at most ~10 traces.
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

    # phylo_mm output (parsed once here; consumed by the branch overlay
    # below and by the header-summary block lower down).
    pm = phylo_mm or {"meta": {}, "branches": []}

    # ── 1c. Driver overlay (phylo_mm) ────────────────────────────────────────
    # One annotation per driver (keyed by its carrier set), each in its own
    # palette colour; its MSA columns are colour-matched (5b).
    _all_tip_names = frozenset(
        parse_leaf_meta(t.name)["sample_hap"] for t in tree.get_terminals()
    )

    def _find_branch_endpoints(rec: dict) -> tuple[float, float, float] | None:
        """Resolve a phylo_mm row to (px, cx, py) endpoints.

        The carrier set can match either half of a bipartition in the drawn
        tree. Try a clade whose tips equal the want set; failing that, find
        the clade equal to the complement and use its sibling's endpoints
        (correct x AND y).

        Using the sibling's actual `pos[id(...)]["y"]` matters: a geometric
        centroid of the want tips disagrees with the recursive
        `(first_child_y + last_child_y)/2` rule that `compute_positions`
        uses, especially for ladder-shaped trees.
        """
        want = frozenset(
            parse_leaf_meta(d)["sample_hap"] for d in rec.get("descendants", ())
        )
        if not want:
            return None
        want_complement = _all_tip_names - want

        def _endpoints_for_clade(c, parent):
            if parent is None:
                return None
            return (pos[id(parent)]["x"], pos[id(c)]["x"], pos[id(c)]["y"])

        # Pass 1: direct match. Pass 2: complement match → resolve sibling.
        for clade in tree.root.find_clades(order="level"):
            if clade.is_terminal() or clade is tree.root:
                continue
            tip_names = frozenset(
                parse_leaf_meta(t.name)["sample_hap"]
                for t in clade.get_terminals()
            )
            parent = next(
                (p for p in tree.root.find_clades() if clade in p.clades),
                None,
            )
            if parent is None:
                continue
            if tip_names == want:
                return _endpoints_for_clade(clade, parent)
            if tip_names == want_complement:
                # Find the sibling clade that holds the want tips directly.
                for sibling in parent.clades:
                    if sibling is clade:
                        continue
                    sib_tips = frozenset(
                        parse_leaf_meta(t.name)["sample_hap"]
                        for t in sibling.get_terminals()
                    )
                    if sib_tips == want:
                        return _endpoints_for_clade(sibling, parent)
                # Fall through — want spans multiple siblings (rare on
                # binary trees). Use the matched clade with the centroid y
                # as a best-effort fallback.
                want_ys = [
                    pos[id(t)]["y"] for t in tree.get_terminals()
                    if parse_leaf_meta(t.name)["sample_hap"] in want
                ]
                cy = sum(want_ys) / len(want_ys) if want_ys else pos[id(clade)]["y"]
                return (pos[id(parent)]["x"], pos[id(clade)]["x"], cy)
        return None

    # ── Per-driver carrier columns ────────────────────────────────────────────
    # One column per contributing variant, ordered as the legend; a filled cell
    # means that haplotype carries that driver's allele. SET-NESTING is visible
    # here and nowhere else: if an inner driver's carriers are a subset of an
    # outer driver's, the inner column's filled rows fall inside the outer's.
    if drivers:
        _sh_to_ly = {
            parse_leaf_meta(names[i])["sample_hap"]: ly[i] for i in range(n)
        }
        for j, drv in enumerate(drivers):
            cy = []
            for d in drv["rec"].get("descendants", ()):
                sh = parse_leaf_meta(d)["sample_hap"]
                if sh in _sh_to_ly:
                    cy.append(_sh_to_ly[sh])
            if not cy:
                continue
            fig.add_trace(go.Bar(
                x=[0.86] * len(cy), base=[j - 0.43] * len(cy), y=cy,
                orientation="h", width=0.9,
                marker=dict(color=drv["color"], line=dict(width=0)),
                hoverinfo="skip", showlegend=False,
                xaxis="x8", yaxis="y1",
            ))

    # ── Driver SHAPE column (right of MSA) + per-driver violins ────────────────
    # A distinct marker symbol per driver ties two panels together: the symbol in
    # the shape column marks which haplotypes CARRY that driver (= feed the
    # carrier side of its violin), and the violin to the right plots those same
    # haps' % actuation (LEFT half) against the rest (RIGHT half). Read a row's
    # symbol, find its point in the matching violin. The Δ (mean_in − mean_out)
    # and joint β are annotated above each violin.
    def _hex_rgba(h: str, a: float) -> str:
        h = h.lstrip("#")
        if len(h) == 6:
            r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
            return f"rgba({r},{g},{b},{a})"
        return h

    vio_ticktext: list[str] = []
    if drivers:
        _sh_ly = {
            parse_leaf_meta(names[i])["sample_hap"]: ly[i]
            for i in range(n) if names[i] not in construct_names
        }
        _sh_act = {
            parse_leaf_meta(names[i])["sample_hap"]: acts[i]
            for i in range(n)
            if names[i] not in construct_names and acts[i] is not None
        }
        all_sh = set(_sh_act)
        driver_carriers = [
            {parse_leaf_meta(d)["sample_hap"] for d in drv["rec"].get("descendants", ())}
            for drv in drivers
        ]
        for j, drv in enumerate(drivers):
            sym, col = drv["symbol"], drv["color"]
            carrier_sh = driver_carriers[j]
            # Comparison group = the LOCAL contrast, not the whole tree. If this
            # driver is NESTED (its carriers ⊂ another driver's), compare against
            # the smallest ENCLOSING driver's OTHER carriers (the sister context);
            # otherwise against everything else. This matches the joint β — a
            # nested reversing variant (opposite-sign β) then reads correctly in
            # the violin instead of being washed out by the distant high/low clade.
            parents = [k for k in range(len(driver_carriers))
                       if k != j and driver_carriers[k] > carrier_sh]
            if parents:
                comp_sh = driver_carriers[
                    min(parents, key=lambda k: len(driver_carriers[k]))
                ] - carrier_sh
            else:
                comp_sh = all_sh - carrier_sh
            # (a) shape column — carriers get the driver's filled symbol; the
            #     local comparison group (what Δ/β is measured against) gets a
            #     small grey open circle.
            car_ly  = [_sh_ly[s] for s in carrier_sh if s in _sh_ly]
            rest_ly = [_sh_ly[s] for s in comp_sh if s in _sh_ly]
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
            # (b) violin — HORIZONTAL, one driver per row (y = driver, stacked),
            #     % actuation on x. Carriers below the row line, comparison above.
            car_x  = [_sh_act[s] * 100 for s in carrier_sh if s in _sh_act]
            rest_x = [_sh_act[s] * 100 for s in comp_sh if s in _sh_act]
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
            _ls = drv["links"]
            vio_ticktext.append(_site_label(min(_ls)) if _ls else f"d{j}")
            _b = drv["rec"].get("beta")
            _bits = []
            # Δ is the LOCAL contrast (carriers − comparison group), so it agrees
            # with the joint β; for a nested reversing driver both come out
            # positive even though the marginal-vs-whole-tree Δ would be ~0.
            if car_x and rest_x:
                _bits.append(f"Δ={(sum(car_x) / len(car_x) - sum(rest_x) / len(rest_x)):+.0f}%")
            if isinstance(_b, float):
                _bits.append(f"β={_b:+.2f}")
            if _bits:
                fig.add_annotation(
                    x=111, y=j, xref="x10", yref="y3", text="  ".join(_bits),
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

    # Helper for compact p-value formatting (None → "NA").
    def _fmt_p(p):
        return f"{p:.2g}" if isinstance(p, float) else "NA"
    def _fmt_beta(b):
        return f"{b:+.2f}" if isinstance(b, float) else "NA"
    def _marginal_delta(rec: dict) -> float | None:
        """Marginal Δ = mean_in − mean_out: carriers vs all other haplotypes."""
        mi, mo = rec.get("mean_in"), rec.get("mean_out")
        if isinstance(mi, float) and isinstance(mo, float):
            return mi - mo
        return None

    # Draw one branch per driver.
    for _gi, drv in enumerate(drivers):
        color = drv["color"]
        rep = drv["rec"]

        # Driver — thick horizontal branch line + annotation. (No
        # vertical bracket: drawing one at the MRCA's x spanning the clade's
        # full tip y-range put a vertical line near the root cutting across
        # many rows, mis-reading as if it were a branch. The thick horizontal
        # segment + the coloured leader lines already mark the clade.)
        seg = _find_branch_endpoints(rep)
        if seg is not None:
            px, cx, py = seg
            links = drv["links"]
            sorted_sites = (
                [_site_label(s) for s in sorted(links, key=int)] if links else []
            )
            linked_str_full = ", ".join(sorted_sites) if sorted_sites else "none"
            if len(sorted_sites) <= 4:
                linked_str_short = (
                    ", ".join(sorted_sites) if sorted_sites else "no SNP"
                )
            else:
                linked_str_short = (
                    ", ".join(sorted_sites[:4])
                    + f", +{len(sorted_sites) - 4} more"
                )
            beta = rep.get("beta")
            ci_lo = rep.get("ci_lower"); ci_hi = rep.get("ci_upper")
            ci_str = (
                f" [{ci_lo:+.2f}, {ci_hi:+.2f}]"
                if isinstance(ci_lo, float) and isinstance(ci_hi, float)
                else ""
            )
            delta = _marginal_delta(rep)
            mi, mo = rep.get("mean_in"), rep.get("mean_out")
            means_str = (
                f"actuation in: {mi:.2f}  vs out: {mo:.2f}<br>"
                if isinstance(mi, float) and isinstance(mo, float)
                else ""
            )
            # Headline = the JOINT (conditional) contribution β and its p from
            # the forward-selection regression: "this variant changes actuation
            # by β, after accounting for the others." A nested variant that
            # switches actuation back off gets a negative β. The marginal
            # carriers-vs-rest Δ is shown alongside in hover.
            head_p = rep.get("p")
            fig.add_trace(go.Scatter(
                x=[px, cx], y=[py, py],
                mode="lines",
                line=dict(color=color, width=3.0),
                hoverinfo="text",
                hovertext=(
                    f"<b>Variant at {linked_str_full}</b>  n_carriers={rep['n_below']}<br>"
                    f"{means_str}"
                    f"<b>β</b> (joint contribution){ci_str} = {_fmt_beta(beta)}<br>"
                    f"<b>p</b> (joint) = {_fmt_p(head_p)}<br>"
                    f"<i>marginal Δ (carriers − rest) = {_fmt_beta(delta)}, "
                    f"BH p_fdr = {_fmt_p(rep.get('p_fdr'))}</i>"
                ),
                name=(
                    f"β={_fmt_beta(beta)} p={_fmt_p(head_p)} "
                    f"pos: {linked_str_short}"
                ),
                legendgroup=f"sig_group_{_gi}",
                # Legend entry comes uniformly from §5c (one per contributing
                # variant) so non-clade variants — whose carriers don't form a
                # branch and so never reach this block — still appear.
                showlegend=False,
                xaxis="x1", yaxis="y1",
            ))
            # NB: the β/p label is NOT drawn on the tree branch — a variant's
            # carriers are often NOT a single clade (homoplasy / recurrent
            # allele), so a tree segment can't represent the strongest variant.
            # The label is anchored to the variant's MSA column instead (§5b
            # below), which always exists. The coloured branch segment here +
            # the colour-matched leader lines just show WHICH haplotypes carry
            # it when they do form a clade.

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

    # ── 2. Leaf labels — coloured by paralog when available ──────────────────
    # Each label is drawn separately so its colour can vary; this is a single
    # Scatter trace with a per-point textfont colour array. Suppressed when
    # row height is too small to fit readable text — identity is then
    # available via hover on the PCLAI / actuation / MSA columns.
    if show_labels:
        # Display the bare sample-hap (last "__"-segment) — paralog/highlight
        # identity is already conveyed by the colour + paralog stripe + chip
        # legend, so the tag prefix would just clutter the label.
        display_labels = [
            construct_label.get(nm) or m["sample_hap"]
            for nm, m in zip(names, metas)
        ]
        fig.add_trace(go.Scatter(
            x=[0.0] * n, y=ly,
            mode="text",
            text=display_labels,
            textposition="middle right",
            textfont=dict(size=label_fs, color=leaf_text_colors),
            hovertext=act_hover,
            hoverinfo="text",
            showlegend=False,
            xaxis="x5", yaxis="y1",
        ))

    # ── 3a. Paralog stripe ────────────────────────────────────────────────────
    # Solid colour bar per leaf in the paralog palette so paralog identity
    # stays readable when row height is too small for label text. Only added
    # when paralogs are actually present — a zero-width axis with a visible
    # full-width Bar trace causes Plotly to fall back and paint the entire
    # figure with the bar colour.
    if has_paralogs:
        paralog_hover = [
            f"<b>{nm}</b><br>paralog: {m['paralog'] or '—'}"
            for nm, m in zip(names, metas)
        ]
        fig.add_trace(go.Bar(
            x=[1] * n, y=ly,
            orientation="h",
            marker=dict(color=leaf_text_colors, line=dict(color="#aaaaaa", width=0.3)),
            width=1.0,
            hovertext=paralog_hover,
            hoverinfo="text",
            showlegend=False,
            xaxis="x7", yaxis="y1",
        ))

    # ── 3b. PCLAI ancestry stripe ─────────────────────────────────────────────
    # Same gate as paralog stripe: only emit the trace when PCLAI hit at least
    # one leaf. Otherwise the bars would render in PCLAI_MISSING grey across
    # the entire figure via Plotly's zero-width-domain fallback.
    if pclai_present:
        fig.add_trace(go.Bar(
            x=[1] * n, y=ly,
            orientation="h",
            marker=dict(color=pclai_rgbs, line=dict(color="#aaaaaa", width=0.3)),
            width=1.0,
            hovertext=pclai_hover,
            hoverinfo="text",
            showlegend=False,
            xaxis="x6", yaxis="y1",
        ))

    # ── 3b. Actuation bars ────────────────────────────────────────────────────
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
    for _ri, (color, label) in enumerate(zip(ACT_COLORS, ACT_LABELS)):
        fig.add_trace(go.Bar(
            x=[None], y=[None],
            orientation="h",
            marker=dict(color=color, line=dict(color="#aaaaaa", width=0.4)),
            name=label,
            showlegend=True, legendrank=100 + _ri,
        ))
    fig.add_trace(go.Bar(
        x=[None], y=[None],
        orientation="h",
        marker=dict(color=ACT_MISSING, line=dict(color="#aaaaaa", width=0.4)),
        name="FIRE missing",
        showlegend=True, legendrank=120,
    ))

    # ── 4b. MSA cell legend chips ─────────────────────────────────────────────
    # The two MSA greys (matches-consensus and gap/N) are easy to confuse, and
    # easy to mistake for the FIRE-missing chip above them. Explicit labels
    # remove ambiguity.
    for color, label in (
        (MSA_INVARIANT, "MSA consensus"),
        (MSA_MISSING,   "MSA gap / N"),
    ):
        fig.add_trace(go.Bar(
            x=[None], y=[None],
            orientation="h",
            marker=dict(color=color, line=dict(color="#aaaaaa", width=0.4)),
            name=label,
            legendgroup="msa",
            showlegend=True,
        ))

    # ── 5. MSA heatmap ────────────────────────────────────────────────────────
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

    # ── 5a. Regulatory-element band markers ───────────────────────────────────
    # Two vertical lines bracketing the RE columns, so it's obvious which part
    # of the (flanked) alignment is the peak itself vs. added assembly context.
    # Skipped when the view is already clipped to the RE (re_only_msa) — the
    # whole panel is the RE, so the brackets would just sit on the edges.
    if re_band is not None and not re_only_msa:
        re_lo, re_hi = re_band
        fig.add_trace(go.Scatter(
            x=[re_lo - 0.5, re_lo - 0.5, None, re_hi + 0.5, re_hi + 0.5],
            y=[Y_RANGE[1], Y_RANGE[0], None, Y_RANGE[1], Y_RANGE[0]],
            mode="lines",
            line=dict(color="#111111", width=1.2, dash="dash"),
            hoverinfo="text",
            hovertext=f"regulatory element: cols {re_lo+1}–{re_hi+1}",
            name=f"RE (cols {re_lo+1}–{re_hi+1})",
            legendgroup="re",
            showlegend=True,
            xaxis="x3", yaxis="y1",
        ))

    # ── 5b. Driver columns (linked_sites) in the MSA, colour-matched ─────────
    # A column linked to several drivers takes the first (lowest p_fdr) colour.
    # Contiguous runs (e.g. a 50-column indel) draw as one outlined box instead
    # of N parallel dashes; singleton columns draw as one dash.
    col_to_color: dict[int, str] = {}
    if drivers:
        for drv in drivers:
            color = drv["color"]
            for s in drv["links"]:
                try:
                    c = int(s) - 1
                except (TypeError, ValueError):
                    continue
                if c not in col_to_color:
                    col_to_color[c] = color
    if col_to_color:
        cols_by_color: dict[str, list[int]] = {}
        for c, color in col_to_color.items():
            cols_by_color.setdefault(color, []).append(c)
        for color, cols in cols_by_color.items():
            cols_sorted = sorted(cols)
            # Split into contiguous runs.
            runs: list[tuple[int, int]] = []
            run_start = cols_sorted[0]
            run_end = cols_sorted[0]
            for c0 in cols_sorted[1:]:
                if c0 == run_end + 1:
                    run_end = c0
                else:
                    runs.append((run_start, run_end))
                    run_start = c0
                    run_end = c0
            runs.append((run_start, run_end))
            link_xs: list = []
            link_ys: list = []
            for r_lo, r_hi in runs:
                if r_lo == r_hi:
                    # Singleton: one vertical dashed line.
                    link_xs += [r_lo, r_lo, None]
                    link_ys += [Y_RANGE[1], Y_RANGE[0], None]
                else:
                    # Multi-column run: dashed rectangle outline.
                    x_lo, x_hi = r_lo - 0.5, r_hi + 0.5
                    link_xs += [x_lo, x_hi, x_hi, x_lo, x_lo, None]
                    link_ys += [Y_RANGE[1], Y_RANGE[1], Y_RANGE[0], Y_RANGE[0], Y_RANGE[1], None]
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

    # ── 5c. Per-variant β/p labels at the TOP of each MSA column ──────────────
    # The headline for every contributing variant lives here (not on the tree):
    # a variant's carriers may not be monophyletic, so the column — which always
    # exists — is the reliable anchor. One label per contributing variant, in
    # its colour, at its representative (left-most) linked column, staggered
    # vertically so neighbouring columns don't collide.
    if drivers:
        _msa_top = Y_RANGE[1]   # data-y of the MSA top (reversed axis)
        for _i, drv in enumerate(drivers):
            rep = drv["rec"]
            positions = sorted(drv["links"], key=int)
            if not positions:
                continue
            c0 = int(positions[0]) - 1   # 0-based column → x3 coordinate
            _beta = _fmt_beta(rep.get("beta"))
            _p = _fmt_p(rep.get("p"))
            # One short line only (β/p) so the box is ~1 row tall and the
            # vertical stagger clears it; the position list lives in the legend
            # and the driver column, so it isn't repeated here.
            fig.add_annotation(
                x=c0, y=_msa_top,
                xref="x3", yref="y1",
                yanchor="bottom",
                yshift=4 + 14 * (_i % 4),
                text=f"<b>β={_beta}</b> p={_p}",
                showarrow=True, arrowhead=0, arrowwidth=0.8,
                arrowcolor=drv["color"], ax=0, ay=-8,
                font=dict(size=8, color=drv["color"]),
                bgcolor="rgba(255,255,255,0.85)",
                bordercolor=drv["color"], borderwidth=0.5,
            )
            # Uniform legend entry per contributing variant — so a variant whose
            # carriers aren't a clade (no branch drawn) still appears.
            _short = _site_label(min(drv["links"])) if drv["links"] else "driver"
            fig.add_trace(go.Scatter(
                x=[None], y=[None], mode="lines",
                line=dict(color=drv["color"], width=3.0),
                name=f"{_short} β={_beta}",
                legendgroup=f"sig_group_{_i}",
                showlegend=True, legendrank=3000 + _i,
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

    # ── Layout ────────────────────────────────────────────────────────────────
    fig.update_layout(
        title=None,   # title is rendered as an annotation so it stays inside the header
        height=fig_height,
        width=fig_width,
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
        xaxis8=dict(
            domain=DRV_X,
            anchor="y1",
            title=dict(text="drivers", font=dict(size=9)) if n_drv else None,
            showticklabels=False, showgrid=False, zeroline=False,
            range=[-0.5, max(0.5, n_drv - 0.5)],
        ),
        xaxis6=dict(
            domain=PCL_X,
            anchor="y1",
            title=dict(text="PCLAI", font=dict(size=9)) if pclai_present else None,
            showticklabels=False, showgrid=False, zeroline=False,
            range=[0, 1.1],
            visible=pclai_present,
        ),
        xaxis7=dict(
            domain=PAR_X,
            anchor="y1",
            title=dict(text="paralog", font=dict(size=9)) if has_paralogs else None,
            showticklabels=False, showgrid=False, zeroline=False,
            range=[0, 1.1],
            visible=has_paralogs,
        ),
        xaxis3=dict(
            domain=MSA_X,
            anchor="y1",
            range=[msa_x0 - 0.5, msa_x1 + 0.5],   # RE-only window unless region mode
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
        # Shape column (carrier markers, leaf-row axis) + violins (% actuation).
        xaxis9=dict(
            domain=MARK_X, anchor="y1",
            range=[-0.5, max(0.5, n_drv - 0.5)],
            title=dict(text="in violin", font=dict(size=8)) if n_drv else None,
            showticklabels=False, showgrid=False, zeroline=False,
            visible=bool(n_drv),
        ),
        xaxis10=dict(
            domain=VIO_X, anchor="y3", range=[-5, 112],
            title=dict(text="% FIRE", font=dict(size=9)) if n_drv else None,
            tickvals=[0, 50, 100], tickfont=dict(size=8),
            showgrid=True, gridcolor="#eeeeee", zeroline=False,
            visible=bool(n_drv),
        ),
        yaxis3=dict(
            domain=VIO_Y3, anchor="x10",
            range=[max(0.7, n_drv - 0.3), -0.7],   # driver 0 at top, stacked down
            tickvals=list(range(len(vio_ticktext))), ticktext=vio_ticktext,
            tickfont=dict(size=8), showgrid=False, zeroline=False,
            visible=bool(n_drv),
        ),
    )

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

    # ── Header: peak title + summary stats (fills the void above tree/act) ──
    valid_acts  = [v for v in acts if v is not None]
    mean_act    = sum(valid_acts) / len(valid_acts) if valid_acts else None
    tree_length = max_depth
    var_sites   = sum(1 for col in range(aln_len) if any(
        row[col] != consensus[col] and row[col] < 4 for row in raw_z
    ))
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
    # H² headline = alpha.V=1 (half-Cauchy(0,1)); alpha.V=0.25 is shown when it
    # differs by >0.05. alpha.V=1000 stays in the TSV only: on weakly identified
    # peaks it inflates H² toward 1.
    meta = pm.get("meta", {})
    h2_main = meta.get("h2_alpha1_mean")
    h2_low  = meta.get("h2_alpha025_mean")
    h2_lo   = meta.get("h2_alpha1_lower")
    h2_hi   = meta.get("h2_alpha1_upper")
    # Convergence flag (alpha.V=1 posterior ESS < threshold in phylo_mm.R).
    # When the MCMC chain didn't converge, the point estimate is unreliable.
    conv_flag = " [ESS low]" if str(meta.get("h2_converged", "")).upper() == "FALSE" else ""
    if isinstance(h2_main, (int, float)):
        # Show the 95% HPD so a point estimate from a nearly-flat posterior
        # isn't mistaken for real heritability. A wide HPD means the data don't
        # identify H² (distinct from prior-sensitivity below): the point
        # estimate is just the median of an uninformative posterior.
        hpd_str, wide_flag = "", ""
        if isinstance(h2_lo, (int, float)) and isinstance(h2_hi, (int, float)):
            hpd_str = f" [{h2_lo:.2f}–{h2_hi:.2f}]"
            if (h2_hi - h2_lo) > 0.5:
                wide_flag = " ⚠ unidentified"
        # Prior-sensitivity flag (α=1 vs α=0.25 disagree): a separate failure mode.
        prior_str = (
            f" (α=0.25→{h2_low:.2f})"
            if isinstance(h2_low, (int, float)) and (h2_main - h2_low) >= 0.05
            else ""
        )
        # Attribution mode: which heritability test (RE-only first, then the
        # whole region) was stable+significant and so gated the variant
        # attribution. `none` = neither, so actuation isn't a heritable function
        # of sequence here and no drivers are reported (the headline H² shown is
        # the unidentified one).
        mode = str(meta.get("attribution_mode", "")).strip()
        if mode == "none":
            def _h2_str(prefix):
                m = meta.get(f"{prefix}_mean")
                lo, hi = meta.get(f"{prefix}_lower"), meta.get(f"{prefix}_upper")
                if not isinstance(m, (int, float)):
                    return "NA"
                hpd = (f" [{lo:.2f}–{hi:.2f}]"
                       if isinstance(lo, (int, float)) and isinstance(hi, (int, float))
                       else "")
                return f"{m:.2f}{hpd} ⚠"
            summary_bits.append(
                f"<b>H²<sub>p</sub></b> not sequence-driven "
                f"(RE {_h2_str('h2_re')}, region {_h2_str('h2_region')}){conv_flag}"
            )
        else:
            scope = "RE" if mode == "RE" else ("region" if mode == "region" else "")
            lbl = f"H²<sub>p</sub>({scope})" if scope else "H²<sub>p</sub>"
            summary_bits.append(
                f"<b>{lbl}</b>≈{h2_main:.2f}{hpd_str}{prior_str}{wide_flag}{conv_flag}"
            )
    # Headline count = the variants the forward regression flagged as
    # contributing (drawn on the tree).
    summary_bits.append(f"<b>contributing variants</b>={n_drv}")
    summary_text = "  ·  ".join(summary_bits)
    # Title left-aligned, summary stats just under it.
    fig.add_annotation(
        text=f"<b>{peak_id}</b>",
        x=0.0, y=1.0, xref="paper", yref="paper",
        xanchor="left", yanchor="top",
        showarrow=False,
        font=dict(size=14, family="Arial", color="#111111"),
    )
    fig.add_annotation(
        text=summary_text,
        x=0.0, y=1.0 - (24.0 / fig_height),
        xref="paper", yref="paper",
        xanchor="left", yanchor="top",
        showarrow=False,
        font=dict(size=11, family="Arial", color="#444444"),
    )

    # ── Paralog legend (only if leaves carry a parsable paralog field) ──────
    # Right-anchored at the MSA right edge so chips never truncate when there
    # are many paralogs or long names; line-break every 4 chips to keep the
    # header tidy.
    if has_paralogs:
        items = [
            f"<span style='color:{c}'><b>■</b> {p}</span>"
            for p, c in paralog_color_map.items()
        ]
        rows = ["  ".join(items[i : i + 4]) for i in range(0, len(items), 4)]
        chips = "<br>".join(rows)
        fig.add_annotation(
            text=chips,
            x=1.0, y=1.0,
            xref="paper", yref="paper",
            xanchor="right", yanchor="top",
            showarrow=False,
            font=dict(size=11, family="Arial"),
            align="right",
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
