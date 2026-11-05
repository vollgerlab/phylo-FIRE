#!/usr/bin/env python3
"""Build a peak -> reference-coordinate table from the union-peaks BED on stdin.

One row per consensus_peak_id with its CHM13 (cons_*) and hg38 (hg38_*) intervals
and `exists_in_hg38`. hg38 coords are constant across a peak's sample-hap rows
(they're the consensus peak's projection), so the first row per peak wins.

Read straight from the union BED so `project_variants` can get hg38 coords
WITHOUT re-running extract_fastas (which would cascade the whole pipeline).

    bgzip -cd table.bed.gz | build_peak_coords.py --out peak_coords.tsv [--peak-ids ids.txt]
"""
import argparse
import sys
from pathlib import Path

# logical -> header name (resolved by NAME, robust to column reordering).
NEEDED = {
    "peak_id": "consensus_peak_id",
    "cons_chr": "cons_chr", "cons_start": "cons_start", "cons_end": "cons_end",
    "hg38_chr": "hg38_chr", "hg38_start": "hg38_start", "hg38_end": "hg38_end",
    "exists": "exists_in_hg38",
}
OUT_COLS = ["peak_id", "cons_chr", "cons_start", "cons_end",
            "hg38_chr", "hg38_start", "hg38_end", "exists_in_hg38"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--peak-ids", type=Path, default=None,
                    help="optional allowlist; only emit these consensus_peak_ids")
    args = ap.parse_args()

    allow = None
    if args.peak_ids:
        allow = {ln.strip() for ln in args.peak_ids.read_text().splitlines() if ln.strip()}

    idx: dict[str, int] | None = None
    seen: dict[str, list[str]] = {}
    for line in sys.stdin:
        line = line.rstrip("\n")
        if idx is None:
            if not line.startswith("#"):
                sys.exit("build_peak_coords: expected a '#'-prefixed header row")
            names = [h.lstrip("#") for h in line.split("\t")]
            colmap = {n: i for i, n in enumerate(names)}
            missing = [v for v in NEEDED.values() if v not in colmap]
            if missing:
                sys.exit(f"build_peak_coords: missing columns {missing}")
            idx = {k: colmap[v] for k, v in NEEDED.items()}
            continue
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) <= max(idx.values()):
            continue
        pid = f[idx["peak_id"]]
        if allow is not None and pid not in allow:
            continue
        if pid in seen:
            continue                    # first row per peak wins (hg38 is constant)
        # sanitize peak id to match on-disk filenames (extract_fastas.sanitize)
        safe = pid.replace("#", "_")
        seen[pid] = [
            safe,
            f[idx["cons_chr"]], f[idx["cons_start"]], f[idx["cons_end"]],
            f[idx["hg38_chr"]], f[idx["hg38_start"]], f[idx["hg38_end"]],
            f[idx["exists"]],
        ]

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as fh:
        fh.write("\t".join(OUT_COLS) + "\n")
        for row in seen.values():
            fh.write("\t".join(row) + "\n")
    print(f"build_peak_coords: wrote {len(seen)} peaks -> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
