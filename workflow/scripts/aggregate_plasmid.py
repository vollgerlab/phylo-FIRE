#!/usr/bin/env python3
"""Concatenate per-peak plasmid-construct TSVs into the single panel deliverable.

One table, `plasmid_constructs.tsv` — every construct for every single-driver
site, with its clonable `re_seq`, the driver + rep support/actuation, and an
`in_random_sample` flag marking the seeded ~N-site subset to actually order. Plus
`plasmid_constructs.fa` — the unique construct sequences (headers carry the
peak/driver/role + in_sample), ready for synthesis.

The per-peak inputs are temp() intermediates; this is the only persisted output.
"""
import argparse
import random
from pathlib import Path

CCOLS = ["peak_id", "driver_id", "driver_beta", "driver_re_cols", "role", "background",
         "carries_driver", "is_synthetic", "collapsed", "source_rep", "source_leaf",
         "source_n_haps", "source_act_mean", "edit", "re_seq", "re_len"]
OUT_COLS = CCOLS + ["in_random_sample"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", required=True, type=Path)
    ap.add_argument("--peaks", required=True, type=Path, help="run manifest")
    ap.add_argument("--out-constructs", required=True, type=Path)
    ap.add_argument("--out-fa", required=True, type=Path)
    ap.add_argument("--n", type=int, default=1000, help="target sites for the random subset")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    peaks = [p.strip() for p in args.peaks.read_text().splitlines() if p.strip()]
    rows = []
    for peak in peaks:
        f = args.indir / f"{peak}.tsv"
        if not f.exists():
            continue
        for ln in f.read_text().splitlines()[1:]:
            if ln.strip():
                rows.append(dict(zip(CCOLS, ln.split("\t"))))

    # seeded random sample of SITES = (peak, driver); flag every construct of a
    # chosen site so `in_random_sample==TRUE` filters to the panel to order.
    sites = sorted({(r["peak_id"], r["driver_id"]) for r in rows})
    rng = random.Random(args.seed)
    chosen = set(rng.sample(sites, min(args.n, len(sites)))) if sites else set()
    for r in rows:
        r["in_random_sample"] = "TRUE" if (r["peak_id"], r["driver_id"]) in chosen else "FALSE"

    args.out_constructs.parent.mkdir(parents=True, exist_ok=True)
    with args.out_constructs.open("w") as fh:
        fh.write("\t".join(OUT_COLS) + "\n")
        for r in rows:
            fh.write("\t".join(r[c] for c in OUT_COLS) + "\n")

    # FASTA of the unique construct sequences (skip collapsed duplicates)
    n_fa = 0
    with args.out_fa.open("w") as fh:
        for r in rows:
            if r["collapsed"] == "TRUE" or not r["re_seq"]:
                continue
            hdr = (f"{r['peak_id']}__{r['driver_id']}__{r['role']} "
                   f"in_sample={r['in_random_sample']} beta={r['driver_beta']} "
                   f"n_haps={r['source_n_haps']} act={r['source_act_mean']}")
            fh.write(f">{hdr}\n{r['re_seq']}\n")
            n_fa += 1

    n_sample = sum(1 for r in rows if r["in_random_sample"] == "TRUE")
    print(f"[plasmid] {len(rows)} constructs over {len(sites)} driver sites "
          f"({len(chosen)} in random-{args.n} sample -> {n_sample} rows); {n_fa} seqs in FASTA")


if __name__ == "__main__":
    main()
