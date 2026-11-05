#!/usr/bin/env python3
"""One-pass FASTA + actuation TSV extraction from union-peaks BED (tab-separated, bgzip'd).

Reads from stdin. Columns are resolved by NAME from the `#`-prefixed header
(robust to reordering across dataset versions), falling back to legacy fixed
indices for a headerless file. Columns used: consensus_peak_id, coverage,
fire_coverage, fire_cov_OV_cov (actuation 0-1), seq, Individual_ID, Haplotype,
is_primary_sample, haplotype_distinct, and the per-sample assembly coordinates
(either explicit asm_chr/asm_start/asm_end, or parsed from
asm_peak_id = '{asm_chr}_{start}_{end}').

Usage:
    bgzip -cd -@ 16 table.bed.gz | python extract_fastas.py --outdir fastas/ [opts]
"""

import argparse
import sys
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path

# Columns are resolved by header name (the union-peaks BED has a `#chrom ...` header).
COLUMN_NAMES = {                       # logical -> header column name
    "peak_id": "consensus_peak_id",
    "coverage": "coverage",
    "fire_coverage": "fire_coverage",
    "actuation": "fire_cov_OV_cov",
    "seq": "seq",
    "individual": "Individual_ID",
    "haplotype": "Haplotype",
    "is_primary": "is_primary_sample",
    "hap_distinct": "haplotype_distinct",
}

# Fallback fixed indices for a headerless legacy file (old cons-coord schema).
LEGACY_IDX = {
    "peak_id": 4, "coverage": 6, "fire_coverage": 7, "actuation": 9,
    "seq": 23, "individual": 25, "haplotype": 26,
    "is_primary": 36, "hap_distinct": 37,
    "_asm": ("split", 11, 12, 13),     # asm_chr, asm_start, asm_end
}


def resolve_columns(header_fields: list[str]) -> dict:
    """Map logical column -> index from a `#`-prefixed header row.

    asm coordinates come from explicit asm_chr/asm_start/asm_end columns or from
    splitting `asm_peak_id = '{asm_chr}_{start}_{end}'`. Records the chosen mode
    under the ``_asm`` key.
    """
    names = [h.lstrip("#") for h in header_fields]
    idx = {n: i for i, n in enumerate(names)}
    cols: dict = {}
    for logical, header_name in COLUMN_NAMES.items():
        if header_name not in idx:
            sys.exit(f"extract_fastas: required column {header_name!r} not in header")
        cols[logical] = idx[header_name]
    if {"asm_chr", "asm_start", "asm_end"} <= idx.keys():
        cols["_asm"] = ("split", idx["asm_chr"], idx["asm_start"], idx["asm_end"])
    elif "asm_peak_id" in idx:
        cols["_asm"] = ("peakid", idx["asm_peak_id"])
    else:
        cols["_asm"] = ("none",)
    return cols


def asm_coords(fields: list[str], cols: dict) -> tuple[str, str, str]:
    """Return (asm_chr, asm_start, asm_end) for a row, per the resolved mode."""
    mode = cols["_asm"]
    if mode[0] == "split":
        return fields[mode[1]], fields[mode[2]], fields[mode[3]]
    if mode[0] == "peakid":
        parts = fields[mode[1]].rsplit("_", 2)   # '{chr}_{start}_{end}'
        if len(parts) == 3:
            return parts[0], parts[1], parts[2]
    return "NA", "NA", "NA"


def sanitize(peak_id: str) -> str:
    """Replace filesystem-unsafe characters in peak IDs."""
    return peak_id.replace("#", "_")


# Consolidated per-haplotype columns (the Stage-2 `haplotypes.parquet`): the same
# fields the per-peak .tsv/.fa carry, in one file so Stage 2 needs no `fastas/` dir.
HAP_COLS = ["peak_id", "sample", "seq", "actuation", "coverage",
            "fire_coverage", "asm_chr", "asm_start", "asm_end"]

# (fasta_header, seq, actuation, coverage, fire_coverage, asm_chr, asm_start, asm_end)
Entry = tuple[str, str, str, str, str, str, str, str]


@dataclass
class Counts:
    """Row tallies reported on stderr after the stream is consumed."""

    rows: int = 0
    skipped_allowlist: int = 0


def parse_row(
    fields: list[str],
    cols: dict,
    peak_id: str,
    primary_only: bool,
    highlight: set[tuple[str, str]],
    highlight_tag: str,
) -> tuple[str, Entry] | None:
    """Return (haplotype_distinct, entry) for a row, or None if it is filtered out."""
    if primary_only and fields[cols["is_primary"]] != "TRUE":
        return None
    seq = fields[cols["seq"]]
    if seq == "NA":
        return None
    sample = f"{fields[cols['individual']]}_{fields[cols['haplotype']]}"
    if (peak_id, sample) in highlight:
        sample = f"{highlight_tag}__{sample}"
    entry = (
        sample, seq, fields[cols["actuation"]], fields[cols["coverage"]],
        fields[cols["fire_coverage"]], *asm_coords(fields, cols),
    )
    return fields[cols["hap_distinct"]], entry


def iter_peaks(
    lines: Iterable[str],
    counts: Counts,
    allowed_peaks: frozenset[str] | None = None,
    primary_only: bool = False,
    highlight: set[tuple[str, str]] | frozenset = frozenset(),
    highlight_tag: str = "highlight",
) -> Iterator[tuple[str, list[Entry], set[str]]]:
    """Group consecutive rows by consensus_peak_id.

    Yields (peak_id, entries, haplotype_distinct values) for every run of rows
    sharing a peak id, including runs whose rows were all filtered out (empty
    entries). The first line resolves the schema: a `#`-prefixed header maps
    columns by name; otherwise the legacy fixed indices apply.
    """
    cols: dict | None = None
    min_fields = 0
    current_peak_id: str | None = None
    entries: list[Entry] = []
    hd_values: set[str] = set()

    for raw_line in lines:
        line = raw_line.rstrip("\n")
        if cols is None:
            if line.startswith("#"):
                cols = resolve_columns(line.split("\t"))
                min_fields = max(v for v in cols.values() if isinstance(v, int)) + 1
                continue
            cols = dict(LEGACY_IDX)
            min_fields = 38
        elif line.startswith("#"):
            continue

        fields = line.split("\t")
        if len(fields) < min_fields:
            continue
        counts.rows += 1
        peak_id = fields[cols["peak_id"]]

        # The peak transition runs before row filters, so a peak with no
        # accepted rows still ends the previous peak.
        if peak_id != current_peak_id:
            if current_peak_id is not None:
                yield current_peak_id, entries, hd_values
            current_peak_id, entries, hd_values = peak_id, [], set()

        if allowed_peaks is not None and peak_id not in allowed_peaks:
            counts.skipped_allowlist += 1
            continue
        parsed = parse_row(fields, cols, peak_id, primary_only, highlight, highlight_tag)
        if parsed is None:
            continue
        haplotype_distinct, entry = parsed
        if haplotype_distinct != "NA":
            hd_values.add(haplotype_distinct)
        entries.append(entry)

    if current_peak_id is not None:
        yield current_peak_id, entries, hd_values


class PeakWriter:
    """Write per-peak FASTA + actuation TSV files that pass the peak filters.

    Tracks the written (sanitized) ids and, when ``collect_haplotypes`` is set,
    accumulates every written row column-major for the consolidated parquet.
    """

    def __init__(
        self,
        outdir: Path,
        min_seqs: int,
        haplotype_distinct: str | None = None,
        collect_haplotypes: bool = False,
    ) -> None:
        self.outdir = outdir
        self.min_seqs = min_seqs
        self.haplotype_distinct = haplotype_distinct
        self.written: list[str] = []
        self.hap_cols: dict[str, list] | None = (
            {k: [] for k in HAP_COLS} if collect_haplotypes else None
        )

    def passes(self, entries: list[Entry], hd_values: set[str]) -> bool:
        """Peak-level filters: non-empty, haplotype_distinct match, >= min_seqs unique."""
        if not entries:
            return False
        if self.haplotype_distinct is not None and self.haplotype_distinct not in hd_values:
            return False
        return len({seq for _, seq, *_ in entries}) >= self.min_seqs

    def write(self, peak_id: str, entries: list[Entry], hd_values: set[str]) -> None:
        """Write one peak's .fa and .tsv if it passes the filters."""
        if not self.passes(entries, hd_values):
            return
        safe_id = sanitize(peak_id)
        with (self.outdir / f"{safe_id}.fa").open("w") as fh:
            for header, seq, *_ in entries:
                fh.write(f">{header}\n{seq}\n")
        with (self.outdir / f"{safe_id}.tsv").open("w") as fh:
            fh.write("sample\tactuation\tcoverage\tfire_coverage\tasm_chr\tasm_start\tasm_end\n")
            for header, _, *rest in entries:
                fh.write("\t".join((header, *rest)) + "\n")
        if self.hap_cols is not None:
            for entry in entries:
                for k, v in zip(HAP_COLS, (safe_id, *entry)):
                    self.hap_cols[k].append(v)
        self.written.append(safe_id)

    def write_parquet(self, path: Path) -> None:
        """Write the accumulated haplotype rows (requires ``collect_haplotypes``)."""
        import polars as pl

        assert self.hap_cols is not None
        path.parent.mkdir(parents=True, exist_ok=True)
        (
            pl.DataFrame(self.hap_cols, schema={c: pl.Utf8 for c in HAP_COLS})
            .with_columns(
                pl.col("actuation").cast(pl.Float64, strict=False),
                pl.col("coverage").cast(pl.Float64, strict=False),
                pl.col("fire_coverage").cast(pl.Float64, strict=False),
            )
            .write_parquet(path)
        )
        print(f"wrote {len(self.hap_cols['peak_id'])} haplotype rows to {path}", file=sys.stderr)

    def write_manifest(self) -> None:
        """Write peak_ids.txt, the authoritative list of peaks produced this run.

        Downstream rules read it instead of globbing the directory, so stale
        files from prior runs never enter the build graph.
        """
        manifest = self.outdir / "peak_ids.txt"
        manifest.write_text("\n".join(sorted(self.written)) + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract per-peak FASTA + actuation TSV files from union-peaks BED on stdin.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--outdir", required=True, type=Path,
        help="Output directory for FASTA and TSV files",
    )
    parser.add_argument(
        "--min-seqs", type=int, default=2,
        help="Minimum unique sequences required to write a FASTA",
    )
    parser.add_argument(
        "--primary-only", action="store_true",
        help="Only include primary samples (is_primary_sample == TRUE)",
    )
    parser.add_argument(
        "--peak-ids", type=Path, default=None,
        help="File of consensus_peak_ids to include (one per line); all others skipped",
    )
    parser.add_argument(
        "--haplotype-distinct", type=str, default=None,
        help=(
            "Only process peaks that have at least one row with this haplotype_distinct "
            "value (e.g. 'Haplotype-selective')"
        ),
    )
    parser.add_argument(
        "--highlight-file", type=Path, default=None,
        help=(
            "Optional TSV (peak_id<TAB>sample_hap, one per line). Matching "
            "leaves get their FASTA header rewritten to "
            "'highlight__{sample_hap}', which the plot script's "
            "paralog-aware leaf-name parser picks up — those leaves then "
            "render with the paralog colour palette + paralog stripe."
        ),
    )
    parser.add_argument(
        "--highlight-tag", type=str, default="highlight",
        help="Tag inserted as the first '__' segment of highlighted leaf names",
    )
    parser.add_argument(
        "--haplotypes-parquet", type=Path, default=None,
        help=(
            "Also write one consolidated parquet (peak_id, sample, seq, "
            "actuation, coverage, fire_coverage, asm_chr/start/end) covering "
            "every written haplotype — the Stage-2 source of actuation + coords "
            "+ RE sequence, replacing the per-peak fastas dir."
        ),
    )
    return parser.parse_args()


def read_peak_ids(path: Path) -> frozenset[str]:
    """Read a one-id-per-line allowlist, ignoring blank lines."""
    return frozenset(line.strip() for line in path.read_text().splitlines() if line.strip())


def read_highlight(path: Path) -> set[tuple[str, str]]:
    """Read (peak_id, sample_hap) pairs from a TSV, skipping blank and `#` lines."""
    highlight: set[tuple[str, str]] = set()
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        pk, _, sm = line.partition("\t")
        if pk and sm:
            highlight.add((pk.strip(), sm.strip()))
    return highlight


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    allowed_peaks = read_peak_ids(args.peak_ids) if args.peak_ids is not None else None
    highlight = read_highlight(args.highlight_file) if args.highlight_file is not None else set()

    counts = Counts()
    writer = PeakWriter(
        args.outdir, args.min_seqs, args.haplotype_distinct,
        collect_haplotypes=args.haplotypes_parquet is not None,
    )
    for peak in iter_peaks(
        sys.stdin, counts, allowed_peaks, args.primary_only, highlight, args.highlight_tag
    ):
        writer.write(*peak)

    if args.haplotypes_parquet is not None:
        writer.write_parquet(args.haplotypes_parquet)
    print(
        f"Processed {counts.rows} rows; wrote {len(writer.written)} FASTA files to {args.outdir}",
        file=sys.stderr,
    )
    if counts.skipped_allowlist:
        print(f"  Skipped {counts.skipped_allowlist} rows not in peak-ids allowlist", file=sys.stderr)
    writer.write_manifest()

    # Each requested peak's FASTA is a DECLARED snakemake output, so a peak that
    # produced nothing (all-NA, or fewer than --min-seqs unique) would leave a
    # missing output and fail the run opaquely. Fail here first with the
    # offending ids, so the fix (drop them from the peak-ids file, or lower
    # --min-seqs) is obvious. At --min-seqs=1 only a truly empty peak trips this.
    if allowed_peaks is not None:
        dropped = sorted({sanitize(p) for p in allowed_peaks} - set(writer.written))
        if dropped:
            sys.exit(
                f"extract_fastas: {len(dropped)} requested peak(s) yielded no "
                f"sequences (all-NA, or < --min-seqs={args.min_seqs} unique): "
                f"{', '.join(dropped[:20])}{' …' if len(dropped) > 20 else ''}. "
                "Remove them from the peak-ids file or lower --min-seqs."
            )


if __name__ == "__main__":
    main()
