#!/usr/bin/env python3
"""One-pass FASTA + actuation TSV extraction from union-peaks BED (tab-separated, bgzip'd).

Reads from stdin. Key columns (0-indexed):
  4  consensus_peak_id
  9  fire_cov_OV_cov  (FIRE actuation fraction, 0-1)
  23 seq
  25 Individual_ID
  26 Haplotype
  36 is_primary_sample
  37 haplotype_distinct

Usage:
    bgzip -cd -@ 16 table.bed.gz | python extract_fastas.py --outdir fastas/ [opts]
"""

import argparse
import sys
from pathlib import Path

# Column indices (0-indexed)
COL_PEAK_ID = 4
COL_COVERAGE = 6
COL_FIRE_COVERAGE = 7
COL_ACTUATION = 9
COL_SEQ = 23
COL_INDIVIDUAL_ID = 25
COL_HAPLOTYPE = 26
COL_IS_PRIMARY = 36
COL_HAPLOTYPE_DISTINCT = 37

MIN_COLS = COL_HAPLOTYPE_DISTINCT + 1


def sanitize(peak_id: str) -> str:
    """Replace filesystem-unsafe characters in peak IDs."""
    return peak_id.replace("#", "_")


def flush_peak(
    peak_id: str,
    entries: list[tuple[str, str, str, str, str]],
    peak_hd_values: set[str],
    outdir: Path,
    min_seqs: int,
    haplotype_distinct_filter: str | None,
    written: list[str],
) -> None:
    """Write one FASTA + actuation TSV if the peak passes the filters.

    entries: list of (fasta_header, sequence, actuation_str, coverage_str, fire_coverage_str) tuples
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
        fh.write("sample\tactuation\tcoverage\tfire_coverage\n")
        for header, _, actuation, coverage, fire_coverage in entries:
            fh.write(f"{header}\t{actuation}\t{coverage}\t{fire_coverage}\n")

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
        "--min-seqs", type=int, default=3,
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
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    allowed_peaks: frozenset[str] | None = None
    if args.peak_ids is not None:
        allowed_peaks = frozenset(
            line.strip()
            for line in args.peak_ids.read_text().splitlines()
            if line.strip()
        )

    written: list[str] = []
    current_peak_id: str | None = None
    current_entries: list[tuple[str, str, str, str, str]] = []
    current_hd_values: set[str] = set()
    n_skipped_allowlist = 0
    n_rows = 0

    for raw_line in sys.stdin:
        line = raw_line.rstrip("\n")
        if line.startswith("#"):
            continue

        fields = line.split("\t")
        if len(fields) < MIN_COLS:
            continue

        n_rows += 1
        peak_id = fields[COL_PEAK_ID]

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
                )
            current_peak_id = peak_id
            current_entries = []
            current_hd_values = set()

        # Peak-level allowlist: skip row collection entirely for non-allowed peaks.
        if allowed_peaks is not None and peak_id not in allowed_peaks:
            n_skipped_allowlist += 1
            continue

        # Row-level filters
        if args.primary_only and fields[COL_IS_PRIMARY] != "TRUE":
            continue

        seq = fields[COL_SEQ]
        if seq == "NA":
            continue

        haplotype_distinct = fields[COL_HAPLOTYPE_DISTINCT]
        if haplotype_distinct != "NA":
            current_hd_values.add(haplotype_distinct)

        individual_id = fields[COL_INDIVIDUAL_ID]
        haplotype = fields[COL_HAPLOTYPE]
        actuation = fields[COL_ACTUATION]
        coverage = fields[COL_COVERAGE]
        fire_coverage = fields[COL_FIRE_COVERAGE]
        current_entries.append((f"{individual_id}_{haplotype}", seq, actuation, coverage, fire_coverage))

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
        )

    print(
        f"Processed {n_rows} rows; wrote {len(written)} FASTA files to {args.outdir}",
        file=sys.stderr,
    )
    if n_skipped_allowlist:
        print(f"  Skipped {n_skipped_allowlist} rows not in peak-ids allowlist", file=sys.stderr)


if __name__ == "__main__":
    main()
