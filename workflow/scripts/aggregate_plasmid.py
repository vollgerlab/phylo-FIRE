#!/usr/bin/env python3
"""Concatenate per-peak plasmid-construct TSVs into the panel deliverables:

  * plasmid_constructs.tsv — every construct (the clonable panel).
  * plasmid_sites.tsv      — one row per (peak, driver): construct count, whether
                            it collapsed to 2, the driver rep vs outgroup rep
                            support + actuation, and an `in_random_sample` flag
                            for the seeded ~N-site subset.
"""
import argparse
import random
from pathlib import Path

CCOLS = ["peak_id", "driver_id", "driver_beta", "driver_re_cols", "role", "background",
         "carries_driver", "is_synthetic", "collapsed", "source_rep", "source_leaf",
         "source_n_haps", "source_act_mean", "edit", "re_seq", "re_len"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", required=True, type=Path)
    ap.add_argument("--peaks", required=True, type=Path, help="run manifest")
    ap.add_argument("--out-constructs", required=True, type=Path)
    ap.add_argument("--out-sites", required=True, type=Path)
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

    args.out_constructs.parent.mkdir(parents=True, exist_ok=True)
    with args.out_constructs.open("w") as fh:
        fh.write("\t".join(CCOLS) + "\n")
        for r in rows:
            fh.write("\t".join(r[c] for c in CCOLS) + "\n")

    # one row per (peak, driver) site
    sites: dict[tuple, dict] = {}
    for r in rows:
        key = (r["peak_id"], r["driver_id"])
        s = sites.setdefault(key, {
            "peak_id": r["peak_id"], "driver_id": r["driver_id"],
            "driver_beta": r["driver_beta"], "n_constructs": 0, "collapsed": "FALSE",
            "driver_n_haps": "NA", "driver_act_mean": "NA",
            "outgroup_n_haps": "NA", "outgroup_act_mean": "NA",
        })
        s["n_constructs"] += 1
        if r["collapsed"] == "TRUE":
            s["collapsed"] = "TRUE"
        if r["background"] == "driver":
            s["driver_n_haps"], s["driver_act_mean"] = r["source_n_haps"], r["source_act_mean"]
        else:
            s["outgroup_n_haps"], s["outgroup_act_mean"] = r["source_n_haps"], r["source_act_mean"]
    # n_constructs is over 4 role-rows; report UNIQUE sequences instead
    for key in sites:
        seqs = {r["re_seq"] for r in rows if (r["peak_id"], r["driver_id"]) == key}
        sites[key]["n_unique_seqs"] = len(seqs)

    # seeded random sample of sites
    rng = random.Random(args.seed)
    all_keys = sorted(sites)
    chosen = set(rng.sample(all_keys, min(args.n, len(all_keys)))) if all_keys else set()

    scols = ["peak_id", "driver_id", "driver_beta", "n_unique_seqs", "collapsed",
             "driver_n_haps", "driver_act_mean", "outgroup_n_haps", "outgroup_act_mean",
             "in_random_sample"]
    with args.out_sites.open("w") as fh:
        fh.write("\t".join(scols) + "\n")
        for key in all_keys:
            s = sites[key]
            s["in_random_sample"] = "TRUE" if key in chosen else "FALSE"
            fh.write("\t".join(str(s[c]) for c in scols) + "\n")

    print(f"[plasmid] {len(rows)} constructs over {len(sites)} driver sites "
          f"({len(chosen)} in the random-{args.n} sample)")


if __name__ == "__main__":
    main()
