#!/usr/bin/env python3
"""Single-pass plasmid construct table for a WHOLE run (no snakemake, no DAG).

Loops the driver peaks, builds the 4-construct design for each single-driver
peak, and writes ONE plasmid_constructs.tsv + .fa. Reads only persisted files
(the mafft alignment, actuation TSV, phylo_mm TSV) -- so it's an additive
post-hoc step on a finished run, exactly like aggregate_hg38.

    pixi run python workflow/scripts/build_plasmid_table.py \
        --aligned-dir  temp/<run>/aligned \
        --fastas-dir   results/<run>/fastas \
        --phylo-mm-dir results/<run>/phylo_mm \
        --peaks        results/<run>/tables/plot_drivers.txt \
        --out-constructs results/<run>/plasmid_fiberseq/plasmid_constructs.tsv \
        --out-fa         results/<run>/plasmid_fiberseq/plasmid_constructs.fa

Peaks whose (temp) alignment no longer exists are skipped and reported.
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_plasmid_constructs import build_peak, COLS


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--aligned-dir", required=True, type=Path)
    ap.add_argument("--fastas-dir", required=True, type=Path)
    ap.add_argument("--phylo-mm-dir", required=True, type=Path)
    ap.add_argument("--peaks", required=True, type=Path, help="driver peak list (one id/line)")
    ap.add_argument("--out-constructs", required=True, type=Path)
    ap.add_argument("--out-fa", required=True, type=Path)
    args = ap.parse_args()

    peaks = [p.strip() for p in args.peaks.read_text().splitlines() if p.strip()]
    all_rows, missing = [], []
    for peak in peaks:
        aln = args.aligned_dir / f"{peak}.fa"
        if not aln.exists():
            missing.append(peak)
            continue
        rows, _ = build_peak(
            aln, args.fastas_dir / f"{peak}.tsv", args.phylo_mm_dir / f"{peak}.tsv", peak
        )
        all_rows.extend(rows)

    args.out_constructs.parent.mkdir(parents=True, exist_ok=True)
    with args.out_constructs.open("w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for r in all_rows:
            fh.write("\t".join(str(r[c]) for c in COLS) + "\n")

    n_fa = 0
    with args.out_fa.open("w") as fh:
        for r in all_rows:
            if r["collapsed"] == "TRUE" or r["re_seq"] == "NA":
                continue
            hdr = (f"{r['peak_id']}__{r['driver_id']}__{r['role']} "
                   f"beta={r['driver_beta']} n_haps={r['source_n_haps']} "
                   f"act={r['source_act_mean']}")
            fh.write(f">{hdr}\n{r['re_seq']}\n")
            n_fa += 1

    sites = {(r["peak_id"], r["driver_id"]) for r in all_rows}
    print(f"[plasmid] {len(all_rows)} constructs over {len(sites)} driver sites; "
          f"{n_fa} unique sequences in FASTA; {len(missing)} peaks skipped "
          f"(no alignment)", file=sys.stderr)
    if missing:
        print("  skipped: " + ",".join(missing), file=sys.stderr)


if __name__ == "__main__":
    main()
