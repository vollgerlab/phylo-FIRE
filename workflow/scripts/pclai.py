"""PCLAI ancestry lookup for tree-plot leaves.

Queries the bgzipped + tabix-indexed HPRC R2 PCLAI BED (CHM13 coordinates,
e.g. resources/pclai/hprc_r2.bed.gz) and returns one annotation per leaf: the
window that overlaps the peak most, else the nearest window within
PCLAI_NEAREST_BP. BED col 4 is ``{sample}/h{hap}_...(pc1,pc2)``; col 9 is the
pre-baked itemRgb used for the stripe.
"""

from pathlib import Path

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
