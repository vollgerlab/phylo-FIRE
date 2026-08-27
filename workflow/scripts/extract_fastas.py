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
from pathlib import Path

# Columns are resolved by header name (the union-peaks BED carries a
# `#chrom ...` header). This survives column reordering across dataset
# versions — e.g. fire_cov_OV_cov moved and asm_chr/start/end were replaced by
# a single asm_peak_id between the cons-coord and asm-coord builds.
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

    asm coordinates come from either explicit asm_chr/asm_start/asm_end columns
    (old schema) or by splitting `asm_peak_id = '{asm_chr}_{start}_{end}'`
    (new schema). Records the chosen mode under the ``_asm`` key.
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


def flush_peak(
    peak_id: str,
    entries: list[tuple[str, str, str, str, str, str, str, str]],
    peak_hd_values: set[str],
    outdir: Path,
    min_seqs: int,
    haplotype_distinct_filter: str | None,
    written: list[str],
    hap_cols: dict[str, list] | None = None,
) -> None:
    """Write one FASTA + actuation TSV if the peak passes the filters.

    entries: list of (fasta_header, sequence, actuation, coverage, fire_coverage,
                      asm_chr, asm_start, asm_end) tuples

    When ``hap_cols`` is given, the written peak's rows are also appended to the
    consolidated haplotypes-parquet accumulator (column-major to keep memory low).
    """
    if not entries:
        return
    # Peak-level haplotype_distinct filter: include peaks with at least one
    # row matching the target classification.
    if haplotype_distinct_filter is not None and haplotype_distinct_filter not in peak_hd_values:
        return
    unique_seqs = {seq for _, seq, *_ in entries}
    if len(unique_seqs) < min_seqs:
        return
    safe_id = sanitize(peak_id)

    fa_path = outdir / f"{safe_id}.fa"
    with fa_path.open("w") as fh:
        for header, seq, *_ in entries:
            fh.write(f">{header}\n{seq}\n")

    tsv_path = outdir / f"{safe_id}.tsv"
    with tsv_path.open("w") as fh:
        fh.write("sample\tactuation\tcoverage\tfire_coverage\tasm_chr\tasm_start\tasm_end\n")
        for header, _, actuation, coverage, fire_coverage, asm_chr, asm_start, asm_end in entries:
            fh.write(
                f"{header}\t{actuation}\t{coverage}\t{fire_coverage}\t"
                f"{asm_chr}\t{asm_start}\t{asm_end}\n"
            )
            if hap_cols is not None:
                for k, v in zip(
                    HAP_COLS,
                    (safe_id, header, _, actuation, coverage, fire_coverage,
                     asm_chr, asm_start, asm_end),
                ):
                    hap_cols[k].append(v)

    written.append(safe_id)


def main() -> None:
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
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    allowed_peaks: frozenset[str] | None = None
    if args.peak_ids is not None:
        allowed_peaks = frozenset(
            line.strip()
            for line in args.peak_ids.read_text().splitlines()
            if line.strip()
        )

    highlight: set[tuple[str, str]] = set()
    if args.highlight_file is not None:
        for line in args.highlight_file.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            pk, _, sm = line.partition("\t")
            if pk and sm:
                highlight.add((pk.strip(), sm.strip()))

    written: list[str] = []
    hap_cols: dict[str, list] | None = (
        {k: [] for k in HAP_COLS} if args.haplotypes_parquet else None
    )
    current_peak_id: str | None = None
    current_entries: list[tuple[str, str, str, str, str, str, str, str]] = []
    current_hd_values: set[str] = set()
    n_skipped_allowlist = 0
    n_rows = 0

    cols: dict | None = None
    min_fields = 0

    for raw_line in sys.stdin:
        line = raw_line.rstrip("\n")

        if cols is None:
            # First line resolves the schema: a `#`-prefixed header maps columns
            # by name; otherwise fall back to the legacy fixed indices.
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

        n_rows += 1
        peak_id = fields[cols["peak_id"]]

        # Detect peak transition — always runs before any row-level filtering so
        # we correctly flush even if the peak had no accepted rows.
        if peak_id != current_peak_id:
            if current_peak_id is not None:
                flush_peak(
                    current_peak_id,
                    current_entries,
                    current_hd_values,
                    args.outdir,
                    args.min_seqs,
                    args.haplotype_distinct,
                    written,
                    hap_cols,
                )
            current_peak_id = peak_id
            current_entries = []
            current_hd_values = set()

        # Peak-level allowlist: skip row collection entirely for non-allowed peaks.
        if allowed_peaks is not None and peak_id not in allowed_peaks:
            n_skipped_allowlist += 1
            continue

        # Row-level filters
        if args.primary_only and fields[cols["is_primary"]] != "TRUE":
            continue

        seq = fields[cols["seq"]]
        if seq == "NA":
            continue

        haplotype_distinct = fields[cols["hap_distinct"]]
        if haplotype_distinct != "NA":
            current_hd_values.add(haplotype_distinct)

        individual_id = fields[cols["individual"]]
        haplotype = fields[cols["haplotype"]]
        actuation = fields[cols["actuation"]]
        coverage = fields[cols["coverage"]]
        fire_coverage = fields[cols["fire_coverage"]]
        asm_chr, asm_start, asm_end = asm_coords(fields, cols)
        sample = f"{individual_id}_{haplotype}"
        if (peak_id, sample) in highlight:
            sample = f"{args.highlight_tag}__{sample}"
        current_entries.append((
            sample, seq, actuation, coverage, fire_coverage,
            asm_chr, asm_start, asm_end,
        ))

    # Flush the final peak
    if current_peak_id is not None:
        flush_peak(
            current_peak_id,
            current_entries,
            current_hd_values,
            args.outdir,
            args.min_seqs,
            args.haplotype_distinct,
            written,
            hap_cols,
        )

    if hap_cols is not None:
        import polars as pl

        args.haplotypes_parquet.parent.mkdir(parents=True, exist_ok=True)
        (
            pl.DataFrame(hap_cols, schema={c: pl.Utf8 for c in HAP_COLS})
            .with_columns(
                pl.col("actuation").cast(pl.Float64, strict=False),
                pl.col("coverage").cast(pl.Float64, strict=False),
                pl.col("fire_coverage").cast(pl.Float64, strict=False),
            )
            .write_parquet(args.haplotypes_parquet)
        )
        print(
            f"wrote {len(hap_cols['peak_id'])} haplotype rows to "
            f"{args.haplotypes_parquet}",
            file=sys.stderr,
        )

    print(
        f"Processed {n_rows} rows; wrote {len(written)} FASTA files to {args.outdir}",
        file=sys.stderr,
    )
    if n_skipped_allowlist:
        print(f"  Skipped {n_skipped_allowlist} rows not in peak-ids allowlist", file=sys.stderr)

    # Authoritative manifest of the peaks produced this run. Downstream
    # snakemake rules read this file instead of globbing the directory, so
    # stale files from prior runs never enter the build graph.
    manifest = args.outdir / "peak_ids.txt"
    manifest.write_text("\n".join(sorted(sanitize(p) for p in written)) + "\n")

    # Each requested peak's FASTA is a DECLARED snakemake output, so a peak that
    # produced nothing (all-NA, or fewer than --min-seqs unique) would leave a
    # missing output and fail the run opaquely. Fail here first with the
    # offending ids, so the fix (drop them from the peak-ids file, or lower
    # --min-seqs) is obvious. At --min-seqs=1 only a truly empty peak trips this.
    if allowed_peaks is not None:
        dropped = sorted({sanitize(p) for p in allowed_peaks} - set(written))
        if dropped:
            sys.exit(
                f"extract_fastas: {len(dropped)} requested peak(s) yielded no "
                f"sequences (all-NA, or < --min-seqs={args.min_seqs} unique): "
                f"{', '.join(dropped[:20])}{' …' if len(dropped) > 20 else ''}. "
                "Remove them from the peak-ids file or lower --min-seqs."
            )


if __name__ == "__main__":
    main()
