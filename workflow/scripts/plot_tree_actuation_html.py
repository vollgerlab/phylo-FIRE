#!/usr/bin/env python3
"""Interactive phylogenetic tree + FIRE actuation + MSA visualization via Plotly.

Layout (left to right):
  tree | leaf labels | PCLAI ancestry stripe | FIRE % actuation | MSA
                                              with sequence logo above

Always writes HTML. Also writes PDF via kaleido (already in pixi deps).

The PCLAI stripe uses HPRC R2's pre-baked PCLAI palette (col 9 of the bgzipped
+ tabix-indexed CHM13-coord BED at resources/pclai/hprc_r2_chm13.bed.gz). The
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
        --pclai-bed resources/pclai/hprc_r2_chm13.bed.gz   # optional
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

    Two formats are supported:
      - paralog peaks:  "LOCUS__PARALOG__SAMPLE_HAP"
      - union peaks:    "SAMPLE_HAP"  (no locus/paralog fields)
    """
    parts = name.split("__")
    if len(parts) == 3:
        return {"locus": parts[0], "paralog": parts[1], "sample_hap": parts[2]}
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

def load_actuation(
    path: Path,
) -> tuple[
    dict[str, float | None],
    dict[str, int | None],
    dict[str, tuple[str, int, int]],
]:
    """Return (actuation, coverage, asm_pos) dicts keyed by sample name.

    asm_pos[sample] = (asm_chr, asm_start, asm_end) when the TSV has those
    columns; absent otherwise (older/preextracted TSVs).
    """
    actuation: dict[str, float | None] = {}
    coverage: dict[str, int | None] = {}
    asm_pos: dict[str, tuple[str, int, int]] = {}
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            name = row["sample"]
            val = row.get("actuation", "NA")
            actuation[name] = None if val in ("NA", "") else float(val)
            cov = row.get("coverage", "NA")
            coverage[name] = None if cov in ("NA", "") else int(cov)
            asm_chr = row.get("asm_chr")
            if asm_chr and asm_chr != "NA":
                try:
                    asm_pos[name] = (asm_chr, int(row["asm_start"]), int(row["asm_end"]))
                except (KeyError, ValueError):
                    pass
    return actuation, coverage, asm_pos


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


def _pansn_chrom(leaf_name: str, asm_chr: str) -> str:
    """Normalize an input-table ``asm_chr`` to the PanSN chrom key used by
    HPRC's asm-coord PCLAI BED (``{sample}#{hap}#{contig}``).

    Most non-HG002 rows already arrive PanSN-formatted (``HG02886#1#CMxxxx.1``)
    and pass through unchanged. The HG002 reference uses parental-suffix style
    (``chr12_PATERNAL``); we strip the suffix and prefix with the leaf's own
    ``{sample}#{hap}#``.
    """
    if "#" in asm_chr:
        return asm_chr
    contig = asm_chr
    for suf in ("_PATERNAL", "_MATERNAL"):
        if contig.endswith(suf):
            contig = contig[: -len(suf)]
            break
    sample_hap = _leaf_to_sample_hap(leaf_name)
    sample, hap = sample_hap.rsplit("_", 1)
    return f"{sample}#{hap}#{contig}"


def _row_to_info(f: list[str], distance_bp: int = 0) -> dict:
    """Convert a parsed tabix row into the per-leaf info dict.

    ``distance_bp`` is 0 for rows that directly overlap the peak, and >0 when
    the row was picked by the nearest-window fallback (Pass 1b in
    :func:`query_pclai`).
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


def _match_pclai_rows(
    rows: list[str],
    start: int,
    end: int,
    leaf_names: list[str],
) -> dict[str, dict]:
    """Pick the tabix row with largest peak overlap per sample-hap.

    Rows that don't overlap the peak interval at all are dropped.
    """
    leaves_by_sample_hap: dict[str, list[str]] = {}
    for nm in leaf_names:
        leaves_by_sample_hap.setdefault(_leaf_to_sample_hap(nm), []).append(nm)

    best_by_sample_hap: dict[str, tuple[int, list[str]]] = {}
    for raw in rows:
        f = raw.split("\t")
        try:
            sample, hap_str = f[3].split("/")[:2]
        except (IndexError, ValueError):
            continue
        sample_hap = f"{sample}_{hap_str.lstrip('h')}"
        if sample_hap not in leaves_by_sample_hap:
            continue
        rs, re_ = int(f[1]), int(f[2])
        overlap = max(0, min(re_, end) - max(rs, start))
        if overlap <= 0:
            continue
        prev = best_by_sample_hap.get(sample_hap)
        if prev is None or overlap > prev[0]:
            best_by_sample_hap[sample_hap] = (overlap, f)

    out: dict[str, dict] = {}
    for sample_hap, (_, f) in best_by_sample_hap.items():
        info = _row_to_info(f, distance_bp=0)
        for leaf in leaves_by_sample_hap[sample_hap]:
            out[leaf] = info
    return out


def _match_pclai_nearest(
    rows: list[str],
    start: int,
    end: int,
    leaf_names: list[str],
    max_distance: int,
) -> dict[str, dict]:
    """Pick the row whose interval is *nearest* to the peak per sample-hap.

    Rows beyond ``max_distance`` bp from the peak interval are dropped, and
    rows that overlap (distance 0) win over any non-overlapping row. Returns
    only entries where the chosen row sits within ``max_distance``; the
    ``distance_bp`` field on each info dict records the gap (0 for overlap).
    """
    leaves_by_sample_hap: dict[str, list[str]] = {}
    for nm in leaf_names:
        leaves_by_sample_hap.setdefault(_leaf_to_sample_hap(nm), []).append(nm)

    best_by_sample_hap: dict[str, tuple[int, list[str]]] = {}
    for raw in rows:
        f = raw.split("\t")
        try:
            sample, hap_str = f[3].split("/")[:2]
        except (IndexError, ValueError):
            continue
        sample_hap = f"{sample}_{hap_str.lstrip('h')}"
        if sample_hap not in leaves_by_sample_hap:
            continue
        rs, re_ = int(f[1]), int(f[2])
        # Distance is 0 when the row overlaps the peak; otherwise the gap
        # between the two intervals (positive integer).
        if re_ <= start:
            dist = start - re_
        elif rs >= end:
            dist = rs - end
        else:
            dist = 0
        if dist > max_distance:
            continue
        prev = best_by_sample_hap.get(sample_hap)
        if prev is None or dist < prev[0]:
            best_by_sample_hap[sample_hap] = (dist, f)

    out: dict[str, dict] = {}
    for sample_hap, (dist, f) in best_by_sample_hap.items():
        info = _row_to_info(f, distance_bp=dist)
        for leaf in leaves_by_sample_hap[sample_hap]:
            out[leaf] = info
    return out


def query_pclai(
    pclai_bed: Path,
    peak_id: str,
    leaf_names: list[str],
    asm_pos: dict[str, tuple[str, int, int]] | None = None,
) -> dict[str, dict]:
    """Look up PCLAI annotation per leaf against the unified bed.

    The merged ``hprc_r2.bed.gz`` carries BOTH coordinate spaces — CHM13
    rows keyed by reference chrom (``chr1``), and asm rows keyed by PanSN
    contig (``{sample}#{hap}#{contig}``). They never collide on chrom keys,
    so a tabix query against either chrom-space hits only rows in that space.

    Lookup proceeds in three passes:
      1. **CHM13 overlap.** For union peaks ``peak_id`` parses as
         ``{chrom}_{start}_{end}`` → one tabix query at that interval. For
         paralog groups (``peak_id="pg_64"``), each leaf carries its own
         sub-peak coords as the first ``__``-segment of its name; leaves are
         grouped by sub-peak and one query runs per distinct sub-peak.
      1b. **CHM13 nearest-within-1Mb.** For leaves still unannotated after
         pass 1 (sample-haps whose assemblies have no PCLAI window directly
         overlapping the peak), a wider tabix query at ``peak ± PCLAI_NEAREST_BP``
         picks the *nearest* PCLAI window per sample-hap. This avoids
         "missing" stripes for haplotypes whose closest annotated window
         simply lies in the immediate flanking region. Rows beyond
         ``PCLAI_NEAREST_BP`` are not considered. Inherited annotations
         carry a positive ``distance_bp`` field.
      2. **Asm fallback.** For leaves still unannotated, if ``asm_pos`` is
         supplied (extract_fastas-emitted ``asm_chr / asm_start / asm_end``
         from the per-peak TSV), each leaf gets one tabix query in PanSN
         chrom space. This recovers HG002-only peaks (no CHM13
         representation) and segmental-duplication paralogs whose CHM13
         ortholog windows have no PCLAI coverage.

    Returns ``{leaf_name: {"rgb", "centroid", "pc1", "pc2", "confidence",
    "distance_bp"}}``. Leaves still unannotated after all passes are absent.
    """
    import pysam   # lazy: only required when PCLAI is wired in
    out: dict[str, dict] = {}

    # Sub-peaks that get queried in CHM13 space, retained so Pass 1b can
    # widen the same intervals without re-parsing leaf names.
    chm13_intervals: list[tuple[str, int, int, list[str]]] = []

    tabix = pysam.TabixFile(str(pclai_bed))
    try:
        # ── Pass 1: CHM13 overlap ──
        parsed = _parse_peak_coords(peak_id)
        if parsed is not None:
            chrom, start, end = parsed
            chm13_intervals.append((chrom, start, end, leaf_names))
            try:
                rows = list(tabix.fetch(chrom, start, end))
            except (ValueError, OSError):
                rows = []
            if rows:
                out.update(_match_pclai_rows(rows, start, end, leaf_names))
        else:
            # Paralog mode: group leaves by sub-peak parsed from name[0].
            by_peak: dict[tuple[str, int, int], list[str]] = {}
            for nm in leaf_names:
                head = nm.split("__", 1)[0]
                coords = _parse_peak_coords(head)
                if coords:
                    by_peak.setdefault(coords, []).append(nm)
            for (chrom, start, end), peak_leaves in by_peak.items():
                chm13_intervals.append((chrom, start, end, peak_leaves))
                try:
                    rows = list(tabix.fetch(chrom, start, end))
                except (ValueError, OSError):
                    continue
                if rows:
                    out.update(_match_pclai_rows(rows, start, end, peak_leaves))

        # ── Pass 1b: CHM13 nearest-within-PCLAI_NEAREST_BP ──
        # Widen each Pass 1 interval and pick the closest window for any
        # sample-hap not yet annotated. This handles haplotypes whose PCLAI
        # call set has no window directly over the peak but does have one
        # nearby (common at SDs and assembly-gap edges).
        for chrom, start, end, peak_leaves in chm13_intervals:
            unannot = [nm for nm in peak_leaves if nm not in out]
            if not unannot:
                continue
            wide_start = max(0, start - PCLAI_NEAREST_BP)
            wide_end   = end + PCLAI_NEAREST_BP
            try:
                rows = list(tabix.fetch(chrom, wide_start, wide_end))
            except (ValueError, OSError):
                continue
            if rows:
                out.update(_match_pclai_nearest(
                    rows, start, end, unannot, PCLAI_NEAREST_BP,
                ))

        # ── Pass 2: asm fallback per leaf ──
        if asm_pos:
            for leaf in leaf_names:
                if leaf in out:
                    continue
                pos = asm_pos.get(leaf)
                if pos is None:
                    continue
                asm_chr, a_start, a_end = pos
                pansn_chrom = _pansn_chrom(leaf, asm_chr)
                try:
                    rows = list(tabix.fetch(pansn_chrom, a_start, a_end))
                except (ValueError, OSError):
                    continue
                if rows:
                    out.update(_match_pclai_rows(rows, a_start, a_end, [leaf]))
    finally:
        tabix.close()
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
    pclai_bed: Path | None = None,
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

    # Per-leaf metadata + colour assignment (paralog-aware).
    metas = [parse_leaf_meta(nm) for nm in names]
    paralog_color_map, leaf_text_colors, has_paralogs = assign_paralog_colors(metas)
    n_unique_seqs = len({alignment.get(nm, "") for nm in names if alignment.get(nm)})

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
    fig_width   = 1500
    target_h    = 2200
    row_h       = max(4, min(18, target_h // max(n, 1)))
    label_fs    = max(7, min(12, row_h - 1))
    show_labels = label_fs >= 8 and row_h >= 9
    fig_height  = max(420, n * row_h + 160)

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
    par_frac        = 0.018 if has_paralogs else 0.0
    pcl_frac        = 0.018 if pclai_present else 0.0
    gap             = 0.005
    msa_frac        = 1.0 - (label_frac + par_frac + pcl_frac + act_frac + gap * 6) - 0.20
    tree_frac       = 1.0 - (label_frac + par_frac + pcl_frac + act_frac + msa_frac + gap * 6)
    # Slide order: TREE | LABEL | PAR | PCL | ACT | MSA. PAR/PCL collapse to
    # zero width when their respective data isn't present; their traces stay
    # in the figure to keep trace indices stable.
    cur = 0.0
    TREE_X  = [cur, cur + tree_frac];   cur = TREE_X[1]  + gap
    LABEL_X = [cur, cur + label_frac];  cur = LABEL_X[1] + gap
    PAR_X   = [cur, cur + par_frac];    cur = PAR_X[1]   + gap
    PCL_X   = [cur, cur + pcl_frac];    cur = PCL_X[1]   + gap
    ACT_X   = [cur, cur + act_frac];    cur = ACT_X[1]   + gap * 2
    MSA_X   = [cur, 1.00]

    # Vertical layout: title strip (top, ~40 px) | freq chart (~50 px) | main.
    # The title strip stays a fixed pixel height regardless of figure tallness,
    # so the title + paralog chips never collide with the freq chart bars.
    TITLE_STRIP_FRAC = min(0.10, 44.0 / fig_height)
    LOGO_BAND_FRAC   = min(0.08, 50.0 / fig_height)
    LOGO_Y           = [1.0 - TITLE_STRIP_FRAC - LOGO_BAND_FRAC, 1.0 - TITLE_STRIP_FRAC]
    MAIN_Y           = [0.00, LOGO_Y[0] - 0.005]

    max_depth = max((pos[id(lf)]["x"] for lf in leaves), default=1.0)

    # Shared y-range — explicit to eliminate auto-padding that breaks alignment
    Y_RANGE = [n - 0.5, -0.5]   # reversed: leaf 0 at top

    fig = go.Figure()

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
        fig.add_trace(go.Scatter(
            x=[0.0] * n, y=ly,
            mode="text",
            text=names,
            textposition="middle right",
            textfont=dict(size=label_fs, color=leaf_text_colors),
            hovertext=act_hover,
            hoverinfo="text",
            showlegend=False,
            xaxis="x5", yaxis="y1",
        ))

    # ── 3a. Paralog stripe ────────────────────────────────────────────────────
    # Solid colour bar per leaf in the paralog palette so paralog identity
    # remains readable when row height is too small for label text. Always
    # added (even when has_paralogs=False); the column collapses to zero
    # width when the trace has no useful colour signal.
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
    # Trace always added so indices stay stable; column collapses to zero
    # width when there's no PCLAI data.
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
        name="FIRE missing",
        legendgroup="act",
        showlegend=True,
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
        margin=dict(l=10, r=150, t=20, b=44),
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
    parser.add_argument("--html",      required=True, type=Path,
                        help="Output HTML path")
    parser.add_argument("--pdf",       required=True, type=Path,
                        help="Output PDF path (via kaleido)")
    parser.add_argument("--peak-id",   type=str, default=None)
    parser.add_argument(
        "--pclai-bed", type=Path, default=None,
        help="Optional bgzipped + tabix-indexed HPRC PCLAI asm-coord BED. "
             "When supplied, draws a per-leaf ancestry stripe colored by HPRC's "
             "pre-baked PCLAI palette.",
    )
    args = parser.parse_args()

    tree                = Phylo.read(str(args.tree), "newick")
    # asm_pos is no longer used since the PCLAI lookup is CHM13-coord-based
    # (one tabix query per peak, keyed by peak_id). The loader keeps the
    # third return for forward-compat with any future asm-coord consumers.
    actuation, coverage, _ = load_actuation(args.actuation)
    alignment           = load_alignment(args.alignment)
    peak_id             = args.peak_id or args.tree.stem

    fig = make_figure(
        tree, actuation, coverage, alignment, peak_id,
        pclai_bed=args.pclai_bed,
    )

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

    # PDF via kaleido. Fail loudly: snakemake declares `pdf` as a required
    # output, so a silent PDF failure produces a confusing "missing outputs"
    # error instead of the actual kaleido stack trace.
    fig.write_image(str(args.pdf), format="pdf", scale=2)
    print(f"Saved {args.pdf}", file=sys.stderr)


if __name__ == "__main__":
    main()
