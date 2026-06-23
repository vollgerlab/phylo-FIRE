#!/usr/bin/env python3
"""Indel-coded alignment for IQ-TREE via Simmons & Ochoterena (2000) simple coding.

Produces a partition Nexus + two PHYLIPs:
  {peak_id}_nucleic.phy  -- the DNA alignment in interleaved-style strict PHYLIP
  {peak_id}_indel.phy    -- binary 0/1/? indel matrix (one column per informative gap)
  {peak_id}.nex          -- IQ-TREE partition file (GTR+I+G on DNA, GTR2+ASC on indel)

Simple indel coding rules (Simmons & Ochoterena 2000):
  - Each distinct contiguous INTERNAL gap interval [start, end) is a binary
    character. Internal = the gap touches neither alignment end.
  - For a given character (gap event G) and taxon T:
      * G outside T's observed span (T's leading/trailing region)  -> state ?
      * T has exactly G                                            -> state 1
      * T has a LARGER gap that strictly contains G                -> state ?
      * T has sequence spanning G (no matching/containing gap)     -> state 0
  - Singleton gaps (count == 1) are non-informative and dropped.
  - Terminal (leading/trailing) gaps are treated as MISSING DATA, not coded:
    the per-haplotype sequence is the assembly substring at the projected peak
    coordinates, so a gap touching an alignment end is indistinguishable from
    the extraction window clipping that haplotype. A taxon is scored `?` for
    any character falling in its terminal / unobserved region.
  - Invariant columns (no `0` or no `1` once `?`s are removed) are dropped so
    the matrix is variable-only. This is required for the GTR2+ASC model:
    ascertainment-bias correction conditions on characters being variable, and
    the un-recorded no-indel positions would otherwise inflate branch lengths.

When zero informative indel characters survive the indel partition is stripped
from the Nexus so IQ-TREE runs on the DNA partition only.

Usage:
    python encode_indels.py --aligned <peak>.fa --peak-id <peak> --out-dir <dir>
"""
import argparse
import sys
from pathlib import Path


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name: str | None = None
    parts: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(parts).upper()
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if name is not None:
        seqs[name] = "".join(parts).upper()
    return seqs


def find_gaps(seq: str) -> list[tuple[int, int]]:
    """Return list of (start, end) half-open intervals for each '-' run."""
    gaps: list[tuple[int, int]] = []
    i = 0
    n = len(seq)
    while i < n:
        if seq[i] == "-":
            j = i
            while j < n and seq[j] == "-":
                j += 1
            gaps.append((i, j))
            i = j
        else:
            i += 1
    return gaps


def observed_span(seq: str) -> tuple[int, int]:
    """Half-open [first_base, last_base + 1) of non-gap sequence.

    Returns (0, 0) for an all-gap row (observes nothing), which makes every
    indel character score `?` for that taxon.
    """
    first = last = None
    for i, c in enumerate(seq):
        if c != "-":
            if first is None:
                first = i
            last = i
    return (0, 0) if first is None else (first, last + 1)


def simple_indel_code(seqs: dict[str, str]) -> tuple[list[str], dict[str, str]]:
    """Apply Simmons-Ochoterena simple indel coding with terminal-gap masking.

    Only INTERNAL gaps (flanked by sequence on both sides) register as indel
    characters; a taxon is scored `?` for any character outside its observed
    span (its leading/trailing region). Invariant columns are dropped so the
    matrix is variable-only (required for GTR2+ASC). See module docstring.

    Returns (taxa_order, {taxon: binary_string}).
    """
    taxa = list(seqs.keys())
    if not taxa:
        return taxa, {}
    L = len(next(iter(seqs.values())))
    all_gaps: dict[str, list[tuple[int, int]]] = {t: find_gaps(seqs[t]) for t in taxa}
    spans: dict[str, tuple[int, int]] = {t: observed_span(seqs[t]) for t in taxa}

    # Candidate indel events: INTERNAL gaps (never touch position 0 or L) seen
    # in >=2 taxa. Terminal gaps are excluded — they are length/window
    # artifacts, not indel events.
    gap_counts: dict[tuple[int, int], int] = {}
    for t in taxa:
        for (gs, ge) in all_gaps[t]:
            if gs == 0 or ge == L:
                continue
            gap_counts[(gs, ge)] = gap_counts.get((gs, ge), 0) + 1
    informative = sorted(g for g, c in gap_counts.items() if c >= 2)

    cols: list[list[str]] = []   # one column (over taxa) per kept character
    for (gs, ge) in informative:
        column: list[str] = []
        for t in taxa:
            os_, oe = spans[t]
            if not (os_ <= gs and ge <= oe):
                column.append("?")          # terminal / unobserved region
                continue
            tg = all_gaps[t]
            if (gs, ge) in tg:
                column.append("1")
            elif any(
                ts <= gs and te >= ge and (ts, te) != (gs, ge) for (ts, te) in tg
            ):
                column.append("?")           # larger containing gap -> ambiguous
            else:
                column.append("0")
        # Keep only variable columns (>=1 '0' AND >=1 '1'); invariant columns
        # would break the GTR2+ASC ascertainment correction.
        if "0" in column and "1" in column:
            cols.append(column)

    codes = {t: "".join(col[i] for col in cols) for i, t in enumerate(taxa)}
    return taxa, codes


def write_phylip(path: Path, taxa: list[str], rows: dict[str, str]) -> None:
    if not taxa:
        path.write_text("0 0\n")
        return
    ncol = len(next(iter(rows.values())))
    name_w = max(len(t) for t in taxa) + 2
    with open(path, "w") as fh:
        fh.write(f" {len(taxa)} {ncol}\n")
        for t in taxa:
            fh.write(f"{t.ljust(name_w)}{rows[t]}\n")


def write_partition_nexus(
    nex_path: Path, peak_id: str, n_dna: int, n_indel: int
) -> None:
    lines = ["#nexus", "begin sets;"]
    lines.append(f"\tcharset part1 = {peak_id}_nucleic.phy: 1-{n_dna};")
    if n_indel > 0:
        lines.append(f"\tcharset part2 = {peak_id}_indel.phy: 1-{n_indel};")
        # GTR2+ASC: the indel matrix is variable-only by construction
        # (invariant columns are dropped in simple_indel_code), so ascertainment
        # correction is required to avoid inflated indel-partition branch lengths.
        lines.append("\tcharpartition mine = GTR+I+G:part1, GTR2+ASC:part2;")
    else:
        lines.append("\tcharpartition mine = GTR+I+G:part1;")
    lines.append("end;")
    nex_path.write_text("\n".join(lines) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--aligned", required=True, type=Path)
    ap.add_argument("--peak-id", required=True)
    ap.add_argument("--out-dir", required=True, type=Path)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    seqs = read_fasta(args.aligned)
    if not seqs:
        sys.exit(f"[encode_indels] no sequences in {args.aligned}")

    taxa = list(seqs.keys())
    n_dna = len(next(iter(seqs.values())))
    write_phylip(args.out_dir / f"{args.peak_id}_nucleic.phy", taxa, seqs)

    _, indel_codes = simple_indel_code(seqs)
    n_indel = len(next(iter(indel_codes.values()))) if indel_codes else 0
    if n_indel > 0:
        write_phylip(args.out_dir / f"{args.peak_id}_indel.phy", taxa, indel_codes)

    write_partition_nexus(
        args.out_dir / f"{args.peak_id}.nex", args.peak_id, n_dna, n_indel
    )

    print(
        f"[encode_indels] {args.peak_id}: DNA cols={n_dna}  indel cols={n_indel}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
