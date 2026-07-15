#!/usr/bin/env python3
"""Concatenate per-peak plasmid-construct TSVs into the single panel deliverable.

One table, `plasmid_constructs.tsv` — EVERY construct for EVERY eligible
(single-driver) site, with its clonable `re_seq` and the driver + rep
support/actuation. No downsampling; the per-peak builder already applies the
filters (single distinct driver, min carriers). Plus `plasmid_constructs.fa` —
the unique construct sequences (headers carry peak/driver/role), ready for
synthesis.

The per-peak inputs are temp() intermediates; this is the only persisted output.
"""
import argparse
from pathlib import Path

CCOLS = ["peak_id", "driver_id", "driver_beta", "driver_re_cols", "role", "background",
         "carries_driver", "is_synthetic", "collapsed", "source_rep", "source_leaf",
         "source_n_haps", "source_act_mean", "edit", "re_seq", "re_len"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", required=True, type=Path)
    ap.add_argument("--peaks", required=True, type=Path, help="run manifest")
    ap.add_argument("--out-constructs", required=True, type=Path)
    ap.add_argument("--out-fa", required=True, type=Path)
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

    args.out_constructs.parent.mkdir(parents=True, exist_ok=True)
    with args.out_constructs.open("w") as fh:
        fh.write("\t".join(CCOLS) + "\n")
        for r in rows:
            fh.write("\t".join(r[c] for c in CCOLS) + "\n")

    # FASTA of the unique construct sequences (skip collapsed duplicates)
    n_fa = 0
    with args.out_fa.open("w") as fh:
        for r in rows:
            if r["collapsed"] == "TRUE" or not r["re_seq"]:
                continue
            hdr = (f"{r['peak_id']}__{r['driver_id']}__{r['role']} "
                   f"beta={r['driver_beta']} n_haps={r['source_n_haps']} "
                   f"act={r['source_act_mean']}")
            fh.write(f">{hdr}\n{r['re_seq']}\n")
            n_fa += 1

    sites = {(r["peak_id"], r["driver_id"]) for r in rows}
    print(f"[plasmid] {len(rows)} constructs over {len(sites)} driver sites; "
          f"{n_fa} unique sequences in FASTA")


if __name__ == "__main__":
    main()
